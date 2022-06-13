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
	template<class T, int d>
	class FluidFreeSurface : public Simulator {
		Typedef_VectorD(d);
	public:
		//Needs to be initialized
		real air_density;//we assume the density of fluid is 1
		VectorD gravity_acc;
		FaceFieldDv<T, d> velocity;
		BoundaryConditionDirect<Field<bool, d>> fixed_bc;
		BoundaryConditionDirect<FaceField<T, d>> vol_bc;
		BoundaryConditionDirect<FaceFieldDv<T, d>> velocity_bc;
		LevelSet<d, PointIntpLinearClamp, HOST> levelset;
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

		void Init(json& j, ImplicitGeometry<d>& geom, Field<bool, d>& fixed, FaceField<real, d>& vol, FaceField<bool, d>& face_fixed, FaceField<real, d>& initial_velocity) {
			air_density = Json::Value<real>(j, "air_density", 1e-3);
			gravity_acc = Json::Value<VectorD>(j, "gravity_acc", Vector<T, d>::Unit(1) * (-9.8));
			velocity = initial_velocity;
			fixed_bc.Init(fixed, fixed);
			velocity_bc.Init(face_fixed, initial_velocity);
			vol_bc.Init(face_fixed, vol);
			levelset.Init(velocity.grid, geom);

			poisson.Init(velocity.grid);
			MG_precond.Allocate_Poisson(velocity.grid);
			MGPCG.Init(&poisson, &MG_precond, false, -1, 1e-6);
		}

		void Update_Poisson_System(Field<bool, d>& fixed, FaceField<T, d>& vol, Field<T, d>& div) {
			Grid<d> grid = velocity.grid;
			fixed.Init(grid);
			vol.Init(grid);
			div.Init(grid);

			//Rule for fixed: if in fixed_bc or phi>=0, set fixed=true, otherwise false
			//Rule for div: if fixed, set to 0. If is the fluid cell in interface, modify with jump
			//otherwise div remain the same

			//first mark all boundary condition of fixed cells
			fixed_bc.Apply(fixed);
			//then mark all air cells to fixed
			fixed.Exec_Nodes(
				[&](const VectorDi cell) {
					if (levelset.phi(cell) >= 0) fixed(cell) = true;
					if (fixed(cell)) div(cell) = 0;
				}
			);
			vol.Calc_Faces(
				[&](const int axis, const VectorDi face)->T {
					//If both two nb cells are invalie, then vol is 0 (outside the boundary)
					//If at least one in two nb cells are fixed, then vol is 0
					//If one in two nb cells are invalid, then vol is just 1/rho
					//If two nb cells are valid, that depends on whether they're
					//all-fluid, all-air or an interface
					real density = 0;
					VectorDi cell0 = face - VectorDi::Unit(axis), cell1 = face;
					if (!vol.grid.Valid(cell0)) std::swap(cell0, cell1);
					if (!vol.grid.Valid(cell0)) return 0;
					//cell0 must be valid
					if (fixed(cell0)) return 0;
					real phi0 = levelset.phi(cell0);
					if (!vol.grid.Valid(cell1)) {
						if (phi0 >= 0) density = air_density;
						else density = 1.0;
					}
					else {
						//then both cells should be all valid
						if (fixed(cell1)) return 0;
						real phi1 = levelset.phi(cell1);
						if (phi0 >= 0 && phi1 >= 0) density = air_density;//all air
						else if (phi0 < 0 && phi1 < 0) density = 1.0;//all fluid
						else {
							//cell0 and cell1 are in different phases
							//theta means the fraction of the phase of cell0
							//that works fine for both phi0<0 and phi0>=0
							T theta = phi0 / (phi0 - phi1);
							T den0 = 1.0, den1 = air_density;
							if (phi0 >= 0) std::swap(den0, den1);
							density = theta * den0 + (1.0 - theta) * den1;
						}
					}
					return 1.0 / density;
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

			return;

			pressure_dev.Init(temp_field_dev.grid);
			auto [iter, res] = MGPCG.Solve(pressure_dev.Data(), temp_field_dev.Data());
			Info("Solve poisson with {} iters and residual {}", iter, res);

			Info("solved pressure: \n{}", pressure_dev);

			////velocity+=grad(p)
			//Exterior_Derivative(temp_velocity_dev, pressure_dev);
			//temp_velocity_dev *= poisson.vol;

			//velocity += temp_velocity_dev;
			//velocity_bc.Apply(velocity);

			//Info("velocity after projection: \n{}", velocity);

			//Info("After projection max velocity {}", velocity.Max_Abs());

			//div_host.Calc_Nodes(
			//	[&](const VectorDi cell) ->T{
			//		if (fixed_host(cell)) return 0;
			//		else return (cell[1] - 7) * (0.04);
			//	}
			//);
			//pressure_dev = div_host;
			//Info("tentative pressure: \n{}", pressure_dev);
			//
			//poisson.Apply(temp_field_dev.Data(), pressure_dev.Data());
			Info("tentative Ap: \n{}", temp_field_dev);
			Info("rhs: {}", rhs);
			temp_field_dev -= rhs;
			//Info("residual: \n{}", temp_field_dev);
		}
	};
}