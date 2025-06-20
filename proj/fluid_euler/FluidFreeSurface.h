//////////////////////////////////////////////////////////////////////////
// Fluid Euler with Free Surface
// Copyright (c) (2022-), Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once

#include "Multigrid.h"
#include "ConjugateGradient.h"
#include "Advection.h"
#include "Simulator.h"
#include "IOFunc.h"
#include "LevelSet.h"
#include "MarchingCubes.h"
#include "GridEulerFunc.h"
#include "Interpolator.h"

namespace Meso {
	template<class T, int d>
	std::tuple<Vector<int, d>,Vector<int, d>,T, T> Face_Neighbor_Cells_And_Values(const Field<T, d>& F, const int axis, const Vector<int, d> face, const T outside_val) {
		Typedef_VectorD(d);
		Vector<int, d> cell0 = face - VectorDi::Unit(axis), cell1= face;
		T val0, val1;
		if (F.grid.Valid(cell0)) val0 = F(cell0);
		else val0 = outside_val;
		if (F.grid.Valid(cell1)) val1 = F(cell1);
		else val1 = outside_val;
		return std::make_tuple(cell0, cell1, val0, val1);
	}

	template<class T, int d>
	class FluidFreeSurface : public Simulator {
		Typedef_VectorD(d);
	public:
		//Needs to be initialized
		//physical quantities
		real air_density;//we assume the density of fluid is 1
		VectorD gravity_acc;
		FaceFieldDv<T, d> velocity;
		//define the system behavior
		Field<uchar, d> cell_type;
		BoundaryConditionDirect<FaceFieldDv<T, d>> velocity_bc;
		LevelSet<d> levelset;
		//utilities
		MaskedPoissonMapping<T, d> poisson;
		VCycleMultigridIntp<T, d> MG_precond;
		ConjugateGradient<T, d> MGPCG;

		//Temporary variables, don't need to be initialized
		FaceFieldDv<T, d> temp_velocity_dev;
		FieldDv<T, d> temp_phi_dev;
		FieldDv<T, d> temp_field_dev;
		FieldDv<T, d> pressure_dev;
		FaceField<T, d> vol_host;
		Field<T, d> div_host;

		void Init(json& j, const ImplicitManifold<d>& geom, Field<uchar, d> &_cell_type, FaceField<real, d>& initial_velocity) {
			air_density = Json::Value<real>(j, "air_density", 1e-3);
			gravity_acc = MathFunc::V<d>(Json::Value<Vector3>(j, "gravity_acc", Vector3::Unit(1) * (-9.8)));
			velocity = initial_velocity;
			cell_type = _cell_type;
			FaceField<bool, d> face_fixed(cell_type.grid);
			face_fixed.Calc_Faces(
				[&](const int axis, const VectorDi face) {
					auto [cell0, cell1, val0, val1] = Face_Neighbor_Cells_And_Values(cell_type, axis, face, uchar(-1));
					if (val0 == 2 || val1 == 2) return true;
					return false;
				}
			);
			velocity_bc.Init(face_fixed, initial_velocity);
			levelset.Init(velocity.grid, geom);

			poisson.Init(velocity.grid);
			MG_precond.Init_Poisson(poisson, 2);
			MGPCG.Init(&poisson, &MG_precond, false, -1, 1e-6);
		}

		void Update_Poisson_System() {
			//Step 1: decide cell types
			cell_type.Calc_Nodes(
				[&](const VectorDi cell) {
					if (cell_type(cell) == 2) {
						return 2;
					}
					else if (levelset.phi(cell) >= 0) return 1;
					else return 0;
				}
			);

			//Step 2: set div to 0 for air and solid cells
			div_host.Init(cell_type.grid);
			div_host.Exec_Nodes(
				[&](const VectorDi cell) {
					if (cell_type(cell) != 0) div_host(cell) = 0;
				}
			);

			//Step 4: set vol, and modify div additionally on the interface
			vol_host.Init(cell_type.grid);
			vol_host.Calc_Faces(
				[&](const int axis, const VectorDi face)->T {
					auto [cell0,cell1,type0, type1]= Face_Neighbor_Cells_And_Values<uchar,d>(cell_type, axis, face, (uchar)(-1));
					
					//order: invalid, fluid,air,solid
					if (type0 > type1 || type1 == (uchar)(-1)) {
						std::swap(cell0, cell1);
						std::swap(type0, type1);
					}
					
					if (type0 == uchar(-1) && type1 == uchar(-1)) return 0;
					if (type0 == uchar(-1) && type1 == 0) return 1.0;
					if (type0 == uchar(-1) && type1 == 1) return 1.0 / air_density;
					if (type0 == uchar(-1) && type1 == 2) return 0;
					if (type0 == 0 && type1 == 0) return 1.0;
					if (type0 == 0 && type1 == 1) {//interface!
						//todo: modify div
						real phi0 = levelset.phi(cell0), phi1 = levelset.phi(cell1);
						T theta = phi0 / (phi0 - phi1);
						T den0 = 1.0, den1 = air_density;
						real density = theta * den0 + (1.0 - theta) * den1;
						return 1.0 / density;
					}
					if (type0 == 0 && type1 == 2) return 0;
					//type0,type1 in {AIR,SOLID}, the result is 0
					return 0;
				}
			);
		}

		virtual real CFL_Time(const real cfl) {
			real dx = velocity.grid.dx;
			real max_vel = GridEulerFunc::Linf_Norm(velocity);
			return dx * cfl / max_vel;
		}


		virtual void Output(DriverMetaData& metadata) {
			std::string vts_name = fmt::format("vts{:04d}.vts", metadata.current_frame);
			bf::path vtk_path = metadata.base_path / bf::path(vts_name);
			VTKFunc::Write_VTS(velocity, vtk_path.string());
			
			VertexMatrix<T, d> verts; ElementMatrix<d> elements;
			Marching_Cubes<T, d, HOST>(verts, elements, levelset.phi);
			std::string surface_name = fmt::format("surface{:04d}.vtu", metadata.current_frame);
			bf::path surface_path = metadata.base_path / bf::path(surface_name);
			VTKFunc::Output_Mesh_As_VTU<T, d>(verts, elements, surface_path.string());
		}

		void Extrapolation(FaceFieldDv<T, d>& velocity) {
			const auto grid = velocity.grid;
			FaceField<bool, d> face_valid_mask(grid);
			face_valid_mask.Calc_Faces(
				[&](const int axis, const VectorDi face) {
					VectorD pos = grid.Face_Center(axis, face);
					return IntpLinearClamp::Value(levelset.phi, pos) < 0;
				}
			);

			std::array<InterpolatorLinearWithMask<T, d, HOST>, d> intp;
			for (int i = 0; i < d; i++) intp[i].Init_Shallow(face_valid_mask.Face_Reference(i));
			FaceField<T, d> velocity_host = velocity;
			velocity_host.Calc_Faces(
				[&](const int axis, const VectorDi face) {
					VectorD pos = grid.Face_Center(axis, face);
					real face_phi = IntpLinearClamp::Value(levelset.phi, pos);
					if (face_phi > 0) {
						VectorD interface_pos = levelset.Closest_Point(pos);
						return intp[axis].Value1(velocity_host.Face_Reference(axis), interface_pos);
					}
					else return velocity_host(axis, face);
				}
			);

			velocity = velocity_host;
		}

		virtual void Advance(DriverMetaData& metadata) {
			real dt = metadata.dt;

			//Advection of levelset
			temp_phi_dev = levelset.phi;
			Advection<IntpLinearClamp, 2>::Advect(dt, temp_field_dev, temp_phi_dev, velocity);
			levelset.phi = temp_field_dev;
			levelset.Fast_Marching(-1);//will calculate whole field

			//Advection of velocity
			Advection<IntpLinearPadding0, 2>::Advect(dt, temp_velocity_dev, velocity, velocity);
			velocity = temp_velocity_dev;
			
			//Add body forces
			velocity += (gravity_acc * dt);
			velocity_bc.Apply(velocity);

			//projection
			//vel_div=div(velocity)
			ExteriorDerivativePadding0::Apply(temp_field_dev, velocity);

			div_host = temp_field_dev;
			Update_Poisson_System();

			poisson.Init(cell_type, vol_host);
			MG_precond.Update_Poisson(poisson);
			temp_field_dev = div_host;

			pressure_dev.Init(temp_field_dev.grid);
			auto [iter, res] = MGPCG.Solve(pressure_dev.Data(), temp_field_dev.Data());
			Info("Solve poisson with {} iters and residual {}", iter, res);

			//Info("solved pressure: \n{}", pressure_dev);

			//velocity+=grad(p)
			ExteriorDerivativePadding0::Apply(temp_velocity_dev, pressure_dev);
			temp_velocity_dev *= poisson.vol;

			velocity += temp_velocity_dev;
			velocity_bc.Apply(velocity);

			Info("After projection max velocity {}", GridEulerFunc::Linf_Norm(velocity));

			Extrapolation(velocity);
			velocity_bc.Apply(velocity);

			//Info("velocity after projection: \n{}", velocity);
		}
	};
}