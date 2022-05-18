//////////////////////////////////////////////////////////////////////////
// Fluid Euler
// Copyright (c) (2022-), Bo Zhu, Mengdi Wang, Yitong Deng
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once

#include "Simulator.h"
#include "EulerParticles.h"
#include "MELPIO.h"

namespace Meso {
	template<int d>
	class MELP : public Simulator {
		Typedef_VectorD(d); Typedef_MatrixD(d);
	public:
		MELP() : Simulator() {}

		EulerParticles<d> e_particles;

		void Init() {
			e_particles.Update_Searcher();
			e_particles.Update_Local_Frames();
			e_particles.Orient_Normals();
		}

		virtual real CFL_Time(const real cfl) {
			real dx = e_particles.dx;
			real max_vel = ArrayFunc::Largest_Norm<VectorD>(e_particles.uRef());
			return dx * cfl / max_vel;
		}
		virtual void Output(const bf::path base_path, const int frame) {
			std::string vts_name = fmt::format("vtu{:04d}.vtu", frame);
			bf::path vtk_path = base_path / bf::path(vts_name);
			Write_MELP_E<d>(e_particles, vtk_path.string());
		}
		virtual void Advance(const int current_frame, const real current_time, const real dt) {
#pragma omp parallel for
			for (int i = 0; i < e_particles.Size(); i++) {
				e_particles.u(i) += -VectorD::Unit(1) * dt;
			}
#pragma omp parallel for
			for (int i = 0; i < e_particles.Size(); i++) {
				e_particles.x(i) += e_particles.u(i) * dt;
			}
		}
	};
}