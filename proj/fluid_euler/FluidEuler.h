//////////////////////////////////////////////////////////////////////////
// Fluid Euler
// Copyright (c) (2022-), Bo Zhu, Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once

#include "Multigrid.h"
#include "ConjugateGradient.h"
#include "Advection.h"
#include "Simulator.h"

namespace Meso {
	template<int d>
	class FluidEuler : public Simulator {
		Typedef_VectorD(d);
	public:
		FaceFieldDv<real, d> velocity;
		BoundaryConditionDirect<FaceFieldDv<real, d>> psi_N;

		FaceFieldDv<real, d> temp_velocity;
		FieldDv<real, d> pressure;
		FieldDv<real, d> vel_div;
		ConjugateGradient<real> MGPCG;
		VCycleMultigrid<real> MG_precond;
		void Init(void) {

		}
		virtual real CFL_Time(const real cfl) {
			real dx = velocity.grid.dx;
			real max_vel = velocity.Max_Abs();
			return dx * cfl / max_vel;
		}
		virtual void Output(const std::string base_path, const std::string frame_path) {

		}
		virtual void Advance(const int current_frame, const real current_time, const real dt) {
			//advection
			SemiLagrangian::Advect(dt, temp_velocity, velocity, velocity);
			psi_N.Apply(temp_velocity);

			//projection
			//vel_div=div(velocity)
			Exterior_Derivative(vel_div, temp_velocity);
			int iter; real res;
			MGPCG.Solve(pressure.Data(), vel_div.Data(), iter, res);

			//velocity+=grad(p)
			Exterior_Derivative(temp_velocity, pressure);
			velocity += temp_velocity;
			psi_N.Apply(velocity);
		}
	};
}