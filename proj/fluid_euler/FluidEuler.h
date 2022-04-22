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
		PoissonMapping<real, d> poisson;
		VCycleMultigrid<real> MG_precond;
		ConjugateGradient<real> MGPCG;
		void Init(Field<bool, d>& fixed, FaceField<real, d>& vol, FaceField<bool, d>& face_fixed, FaceField<real, d>& initial_velocity) {
			velocity.Deep_Copy(initial_velocity);
			psi_N.Init(face_fixed, initial_velocity);
			temp_velocity.Init(velocity.grid);
			pressure.Init(velocity.grid);
			vel_div.Init(velocity.grid);

			poisson.Init(velocity.grid, vol, fixed);
			MG_precond.Init_Poisson(poisson, 2, 2);
			MGPCG.Init(&poisson, &MG_precond, false, -1, 1e-5);
		}
		virtual real CFL_Time(const real cfl) {
			real dx = velocity.grid.dx;
			real max_vel = velocity.Max_Abs();
			Info("cfl={}, dx={}, max_vel={}, cfl time {}", cfl, dx, max_vel, dx * cfl / max_vel);
			return dx * cfl / max_vel;
		}
		virtual void Output(const std::string base_path, const std::string frame_path) {

		}
		virtual void Advance(const int current_frame, const real current_time, const real dt) {
			FieldDv<real, d> vel0(velocity.grid.Face_Grid(0), velocity.face_data[0]);
			Info("velocity 0 before advection: \n{}\n", vel0);


			//advection
			SemiLagrangian::Advect(dt, temp_velocity, velocity, velocity);
			psi_N.Apply(temp_velocity);

			Info("velocity 0 after advection: \n{}\n", vel0);

			//projection
			//vel_div=div(velocity)
			Exterior_Derivative(vel_div, temp_velocity);

			VectorDi cell0 = VectorFunc::Vi<d>(0, 0, 0);
			for (int axis = 0; axis < d; axis++) {
				VectorDi face0 = cell0, face1 = cell0 + VectorDi::Unit(axis);
				Info("cell {} incident face {},{} velocity {}", cell0, axis, face0, velocity.Get(axis, face0));
				Info("cell {} incident face {},{} velocity {}", cell0, axis, face1, velocity.Get(axis, face1));
			}

			Info("max abs vel_div before projection: {}", vel_div.Max_Abs());
			Info("vel_div before projection: \n{}\n", vel_div);

			int iter; real res;
			MGPCG.Solve(pressure.Data(), vel_div.Data(), iter, res);

			//velocity+=grad(p)
			Exterior_Derivative(temp_velocity, pressure);

			Info("solved pressure: \n{}\n", pressure);

			velocity += temp_velocity;
			psi_N.Apply(velocity);

			Exterior_Derivative(vel_div, velocity);
			Info("max div abs after projection: {}", vel_div.Max_Abs());

			exit(0);
		}
	};
}