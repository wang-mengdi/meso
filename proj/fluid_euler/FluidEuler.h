//////////////////////////////////////////////////////////////////////////
// Fluid Euler
// Copyright (c) (2022-), Bo Zhu, Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once

#include "Multigrid.h"
#include "ConjugateGradient.h"
#include "Advection.h"

namespace Meso {
	template<int d>
	class FluidEuler : public Simulator {
		Typedef_VectorD(d);
	public:
		FaceFieldDv<real, d> velocity;
		BoundaryConditionDirect<FaceFieldDv<real, d>> psi_N;
		FieldDv<real, d> pressure;
		FieldDv<real, d> vel_div;
		FaceFieldDv<real, d> gradp;
		ConjugateGradient<real> MGPCG;
		VCycleMultigrid<real> MG_precond;
		void Init(void) {

		}
		virtual void Output(const std::string base_path, const std::string frame_path) {

		}
		virtual void Advance(const int current_frame, const real current_time, const real dt) {
			//advection
			SemiLagrangian::Advect(dt, velocity, velocity, psi_N);
			//projection
			//vel_div=div(velocity)
			Exterior_Derivative(vel_div, velocity);
			MGPCG.Apply(pressure.Data(), vel_div.Data());
			//velocity+=grad(p)
			Exterior_Derivative(gradp, pressure);
		}
	};
}