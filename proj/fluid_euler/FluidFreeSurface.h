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

namespace Meso {

	enum CellType : int {
		INVALID = -1,
		FLUID,
		AIR,
		SOLID
	};

	template<class T, int d>
	std::tuple<Vector<int, d>, Vector<int, d>, T, T> Face_Neighbor_Cells_And_Values(const Field<T, d>& F, const int axis, const Vector<int, d> face, const T outside_val) {
		Typedef_VectorD(d);
		VectorDi cell0 = face - VectorDi::Unit(axis), cell1 = face;
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
		Field<CellType, d> cell_type;
		BoundaryConditionDirect<FaceFieldDv<T, d>> velocity_bc;
		LevelSet<d, PointIntpLinearClamp, HOST> levelset;
		//utilities
		MaskedPoissonMapping<T, d> poisson;
		VCycleMultigridIntp<T, d> MG_precond;
		ConjugateGradient<T> MGPCG;

		//Temporary variables, don't need to be initialized
		FaceFieldDv<T, d> temp_velocity_dev;
		FieldDv<T, d> temp_phi_dev;
		FieldDv<T, d> temp_field_dev;
		FieldDv<T, d> pressure_dev;
		Field<bool, d> fixed_host;
		FaceField<T, d> vol_host;
		Field<T, d> div_host;

		void Init(json& j, ImplicitGeometry<d>& geom, Field<CellType, d> _cell_type, FaceField<bool, d>& face_fixed, FaceField<real, d>& initial_velocity) {
			air_density = Json::Value<real>(j, "air_density", 1e-3);
			gravity_acc = Json::Value<VectorD>(j, "gravity_acc", Vector<T, d>::Unit(1) * (-9.8));
			velocity = initial_velocity;
			cell_type = _cell_type;
			velocity_bc.Init(face_fixed, initial_velocity);
			levelset.Init(velocity.grid, geom);

			poisson.Init(velocity.grid);
			MG_precond.Allocate_Poisson(velocity.grid);
			MGPCG.Init(&poisson, &MG_precond, false, -1, 1e-6);
		}

		void Update_Poisson_System(Field<bool, d>& fixed, FaceField<T, d>& vol, Field<T, d>& div) {
			//Step 1: decide cell types
			cell_types.Calc_Nodes(
				[&](const VectorDi cell) {
					if (cell_type(cell) == SOLID) {
						return SOLID;
					}
					else if (levelset.phi(cell) >= 0) return AIR;
					else return FLUID;
				}
			);

			//Step 2: decide fixed from cell types
			fixed.Init(velocity.grid);
			fixed.Calc_Nodes(
				[&](const VectorDi cell) {
					if (cell_type(cell) == FLUID) return false;
					else return true;
				}
			);

			//Step 3: set div to 0 for air and solid cells
			div.Init(velocity.grid);
			div.Exec_Nodes(
				[&](const VectorDi cell) {
					if (cell_type(cell) != FLUID) div(cell) = 0;
				}
			);

			//Step 4: set vol, and modify div additionally on the interface
			vol.Init(grid);
			vol.Calc_Faces(
				[&](const int axis, const VectorDi face)->T {
					auto [cell0, cell1, type0, type1] = Face_Neighbor_Cells_And_Values(cell_type, axis, face, INVALID);
					
					//order: invalid, fluid,air,solid
					if (type0 > type1) {
						std::swap(cell0, cell1);
						std::swap(type0, type1);
					}
					
					if (type0 == INVALID && type1 == INVALID) return 0;
					if (type0 == INVALID && type1 == FLUID) return 1.0;
					if (type0 == INVALID && type1 == AIR) return 1.0 / air_density;
					if (type0 == INVALID && type1 == SOLID) return 0;
					if (type0 == FLUID && type1 == FLUID) return 1.0;
					if (type0 == FLUID && type1 == AIR) {//interface!
						//todo: modify div
						real phi0 = levelset.phi(cell0), phi1 = levelset.phi(cell1);
						T theta = phi0 / (phi0 - phi1);
						T den0 = 1.0, den1 = air_density;
						real density = theta * den0 + (1.0 - theta) * den1;
						return 1.0 / density;
					}
					if (type0 == FLUID && type1 == SOLID) return 0;
					//type0,type1 in {AIR,SOLID}, the result is 0
					return 0;
				}
			);
		}

		virtual real CFL_Time(const real cfl) {
			real dx = velocity.grid.dx;
			real max_vel = velocity.Max_Abs();
			return dx * cfl / max_vel;
		}

		virtual void Output(DriverMetaData& metadata) {
			std::string vts_name = fmt::format("vts{:04d}.vts", metadata.current_frame);
			bf::path vtk_path = metadata.base_path / bf::path(vts_name);
			VTKFunc::Write_VTS(velocity, vtk_path.string());
			
			//VertexMatrix<T, d> verts; ElementMatrix<d> elements;
			//Marching_Cubes<T, d, HOST>(verts, elements, levelset.phi);
			//std::string obj_name = fmt::format("surface{:04d}.obj", metadata.current_frame);
			//bf::path obj_path = metadata.base_path / bf::path(obj_name);
			//OBJFunc::Write_OBJ(obj_path.string(), verts, elements);
		}

		virtual void Advance(DriverMetaData& metadata) {
			real dt = metadata.dt;

			//Advection of levelset
			temp_phi_dev = levelset.phi;
			SemiLagrangian<IntpLinearClamp>::Advect(dt, temp_field_dev, temp_phi_dev, velocity);
			levelset.phi = temp_field_dev;
			//levelset.Fast_Marching(-1);//will calculate whole field

			//Advection of velocity
			SemiLagrangian<IntpLinearPadding0>::Advect(dt, temp_velocity_dev, velocity, velocity);
			velocity = temp_velocity_dev;
			
			//Add body forces
			velocity += (gravity_acc * dt);
			velocity_bc.Apply(velocity);

			//projection
			//vel_div=div(velocity)
			Exterior_Derivative(temp_field_dev, velocity);

			div_host = temp_field_dev;
			Update_Poisson_System(fixed_host, vol_host, div_host);

			poisson.Init(fixed_host, vol_host);
			MG_precond.Update_Poisson(poisson, 2, 2);
			temp_field_dev = div_host;

			Info("phi: \n{}", levelset.phi);
			Info("v: \n{}", velocity);
			Info("fixed_host: \n{}", fixed_host);
			Info("vol: \n{}", vol_host);
			Info("rhs: \n{}", temp_field_dev);
			Info("rhs max: {}", temp_field_dev.Max_Abs());
			FieldDv<T, d> rhs = temp_field_dev;
			//Info("saved rhs: \n{}", rhs);
			Info("rhs data ptr: {}", (void*)rhs.Data_Ptr());
			Info("temp_field_dev data ptr: {}", (void*)temp_field_dev.Data_Ptr());

			pressure_dev.Init(temp_field_dev.grid);
			auto [iter, res] = MGPCG.Solve(pressure_dev.Data(), temp_field_dev.Data());
			Info("Solve poisson with {} iters and residual {}", iter, res);

			Info("solved pressure: \n{}", pressure_dev);

			//velocity+=grad(p)
			Exterior_Derivative(temp_velocity_dev, pressure_dev);
			temp_velocity_dev *= poisson.vol;

			velocity += temp_velocity_dev;
			velocity_bc.Apply(velocity);

			Info("velocity after projection: \n{}", velocity);

			Info("After projection max velocity {}", velocity.Max_Abs());

			div_host.Calc_Nodes(
				[&](const VectorDi cell) ->T{
					if (fixed_host(cell)) return 0;
					else return (cell[1] - 7) * (0.04);
				}
			);
			pressure_dev = div_host;
			Info("tentative pressure: \n{}", pressure_dev);
			
			poisson.Apply(temp_field_dev.Data(), pressure_dev.Data());
			Info("tentative Ap: \n{}", temp_field_dev);
			Info("rhs: {}", rhs);
			temp_field_dev -= rhs;
			Info("residual: \n{}", temp_field_dev);
		}
	};
}