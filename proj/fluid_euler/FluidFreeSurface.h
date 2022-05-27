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

namespace Meso {
	template<class T, int d>
	class FluidFreeSurface : public Simulator {
		Typedef_VectorD(d);
	public:
		T air_density;//we assume the density of fluid is 1
		Vector<T, d> gravity_acc;
		FaceFieldDv<T, d> velocity;
		BoundaryConditionDirect<Field<bool, d>> fixed_bc;
		BoundaryConditionDirect<FaceField<T, d>> vol_bc;
		BoundaryConditionDirect<FaceFieldDv<T, d>> velocity_bc;
		LevelSet<d, PointIntpLinearClamp, HOST> levelset;

		FaceFieldDv<T, d> temp_velocity_dev;
		FieldDv<T, d> temp_phi_dev;
		FieldDv<T, d> temp_field_dev;
		FieldDv<T, d> pressure_dev;
		Field<bool, d> fixed_host;
		FaceField<T, d> vol_host;
		Field<T, d> div_host;
		MaskedPoissonMapping<T, d> poisson;
		VCycleMultigridIntp<T, d> MG_precond;
		ConjugateGradient<T> MGPCG;


		void Update_Poisson_System(Field<bool, d>& fixed, FaceField<T, d>& vol, Field<T, d>& div) {
			Grid<d> grid = velocity.grid;
			fixed.Init(grid);
			vol.Init(grid);
			div.Init(grid);

			//first mark all boundary condition of fixed cells
			fixed_bc.Apply(fixed);
			//then mark all air cells to fixed
			fixed.Exec_Nodes(
				[&](const VectorDi cell) {
					if (levelset.phi(cell) >= 0) fixed(cell) = true;
				}
			);
			vol.Exec_Faces(
				[&](const int axis, const VectorDi face) {
					T density = 0;
					VectorDi cell0 = face - VectorDi::Unit(axis), cell1 = face;
					//there must be at least 1 valid cell
					if (!vol.grid.Valid(cell0)) std::swap(cell0, cell1);
					T phi0 = levelset.phi(cell0), phi1 = levelset.phi(cell1);
					if (!vol.grid.Valid(cell1)) {
						if (phi0 >= 0) density = air_density;
						else density = 1.0;
					}
					else {
						//then both cells should be all valid
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
					vol(axis, face) = 1.0 / density;
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
		}

		virtual void Advance(DriverMetaData& metadata) {
			real dt = metadata.dt;

			//Advection of levelset
			temp_phi_dev = levelset.phi;
			SemiLagrangian<IntpLinearClamp>::Advect(dt, temp_field_dev, temp_phi_dev, velocity);
			levelset.phi = temp_field_dev;
			levelset.Fast_Marching(-1);//will calculate whole field

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

			int iter; real res;
			MGPCG.Solve(pressure_dev.Data(), temp_field_dev.Data(), iter, res);
			Info("Solve poisson with {} iters and residual {}", iter, res);

			//velocity+=grad(p)
			Exterior_Derivative(temp_velocity_dev, pressure_dev);
			temp_velocity_dev *= poisson.vol;

			velocity += temp_velocity_dev;
			velocity_bc.Apply(velocity);

			Exterior_Derivative(temp_field_dev, velocity);
			Info("After projection max div {}", temp_field_dev.Max_Abs());
		}
	};
}