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
	class FluidEuler : public Simulator {
		Typedef_VectorD(d);
	public:
		Vector<T, d> gravity_acc;
		FaceFieldDv<T, d> velocity;
		BoundaryConditionDirect<FaceFieldDv<T, d>> psi_N;
		LevelSet<d, PointIntpLinearClamp, HOST> levelset;

		FaceFieldDv<T, d> temp_velocity;
		Field<T, d> temp_phi;
		FieldDv<real, d> vel_div;
		FieldDv<real, d> pressure;
		MaskedPoissonMapping<real, d> poisson;
		VCycleMultigridIntp<real, d> MG_precond;
		ConjugateGradient<real> MGPCG;

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
			SemiLagrangian<IntpLinearClamp>::Advect(dt, temp_phi, levelset.phi, velocity);
			levelset.phi = temp_phi;
			levelset.Fast_Marching(-1);//will calculate whole field

			//Advection of velocity
			SemiLagrangian<IntpLinearPadding0>::Advect(dt, temp_velocity, velocity, velocity);
			velocity = temp_velocity;
			
			//Add body forces
			velocity += (gravity_acc * dt);
			psi_N.Apply(velocity);

			//projection
			//vel_div=div(velocity)
			Exterior_Derivative(vel_div, velocity);

			int iter; real res;
			MGPCG.Solve(pressure.Data(), vel_div.Data(), iter, res);
			Info("Solve poisson with {} iters and residual {}", iter, res);

			//velocity+=grad(p)
			Exterior_Derivative(temp_velocity, pressure);
			temp_velocity *= poisson.vol;

			velocity += temp_velocity;
			psi_N.Apply(velocity);

			Exterior_Derivative(vel_div, velocity);
			Info("After projection max div {}", vel_div.Max_Abs());
		}
	};
}