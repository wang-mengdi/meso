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
		FaceFieldDv<real, d> velocity;
		BoundaryConditionDirect<FaceFieldDv<real, d>> psi_N;

		virtual real CFL_Time(const real cfl) {
			real dx = velocity.grid.dx;
			real max_vel = velocity.Max_Abs();
			return dx * cfl / max_vel;
		}
		virtual void Output(const bf::path base_path, const int frame) {
			std::string vts_name = fmt::format("vts{:04d}.vts", frame);
			bf::path vtk_path = base_path / bf::path(vts_name);
			VTKFunc::Write_VTS(velocity, vtk_path.string());
		}
		virtual void Advance(const int current_frame, const real current_time, const real dt) {
			//Advection
			//SemiLagrangian::Advect
		}
	};
}