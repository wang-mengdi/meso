//////////////////////////////////////////////////////////////////////////
// FluidSPH
// Copyright (c) (2022-), Yitong Deng
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once

#include "Simulator.h"
#include "SPHParticles.h"
#include "SPHIO.h"
#include "SPH.h"

namespace Meso {

	template<int d>
	class FluidSPH : public Simulator {
		Typedef_VectorD(d); Typedef_MatrixD(d);
	public:
		FluidSPH() : Simulator() {}

		SPHParticles<d> particles;

		void Init() {
			particles.Update_Searcher();
			Register_Nbs();
		}

		virtual real CFL_Time(const real cfl) {
			real dx = particles.dx;
			real max_vel = ArrayFunc::Largest_Norm<VectorD>(particles.uRef());
			return dx * cfl / max_vel;
		}
		virtual void Output(const bf::path base_path, const int frame) {
			std::string vts_name = fmt::format("vtp{:04d}.vtp", frame);
			bf::path vtk_path = base_path / bf::path(vts_name);
			Write_SPH<d>(particles, vtk_path.string());
		}

		real Radius(const int i) const {
			return 4. * particles.dx;
		}

		void Register_Nbs(void) {
			if (particles.nbs_info.size() != particles.Size()) particles.nbs_info.resize(particles.Size());
#pragma omp parallel for
			for (int i = 0; i < particles.Size(); i++) {
				Array<int> nbs; particles.nbs_searcher->Find_Nbs(particles.x(i), Radius(i), nbs);
				particles.nbs_info[i].set(nbs);
			}
		}

		virtual void Advance(const int current_frame, const real current_time, const real dt) {
			Info("Advancing");
		}
	};
}