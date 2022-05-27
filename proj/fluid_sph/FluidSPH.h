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
#include "AnalyticalBoundary.h"

namespace Meso {

	template<int d>
	class FluidSPH : public Simulator {
		Typedef_VectorD(d); Typedef_MatrixD(d);
	public:
		FluidSPH() : Simulator() {}

		SPHParticles<d> particles;
		AnalyticalBoundary<d> boundary;

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

		int Boundary(const int i) const {
			return particles.B(i);
		}

		int Has_Boundary(void) const {
			return boundary.Available();
		}
			
		void Register_Nbs(void) {
			if (particles.nbs_info.size() != particles.Size()) particles.nbs_info.resize(particles.Size());
#pragma omp parallel for
			for (int i = 0; i < particles.Size(); i++) {
				Array<int> nbs; particles.nbs_searcher->Find_Nbs(particles.x(i), Radius(i), nbs);
				particles.nbs_info[i].set(nbs);
			}
		}

		void Body_Force(const real dt) {
#pragma omp parallel for
			for (int i = 0; i < particles.Size(); i++) {
				particles.acc(i) += -VectorD::Unit(0);
			}
		}

		void Update_Position(const real dt) {
#pragma omp parallel for
			for (int i = 0; i < particles.Size(); i++) {
				particles.u(i) += dt * particles.acc(i);
			}
#pragma omp parallel for
			for (int i = 0; i < particles.Size(); i++) {
				if (Boundary(i)) particles.u(i) *= 0.;
			}
#pragma omp parallel for
			for (int i = 0; i < particles.Size(); i++) {
				particles.x(i) += dt * particles.u(i);
			}
		}

		void Prepare_Advance(void) {
#pragma omp parallel for
			for (int i = 0; i < particles.Size(); i++) {
				particles.acc(i) *= 0.;
			}
		}

		void Update_Phis(void)
		{
#pragma omp parallel for
			for (int i = 0; i < particles.Size(); i++) {
				real phi; VectorD normal;
				bool flg = boundary.Get_Nearest_Boundary(particles.x(i), phi, normal);
				particles.bnd_phi(i) = phi;
				particles.bnd_n(i) = normal;
			}
		}

		void Enforce_Boundary(const real dt) {
			Update_Phis();
#pragma omp parallel for
			for (int i = 0; i < particles.Size(); i++) {
				if (!Boundary(i)) {
					real phi = particles.bnd_phi(i);
					if (phi > 0) {
						//particles.B(i) = 1;
						VectorD normal = particles.bnd_n(i);
						VectorD normal_velocity = particles.u(i).dot(normal) * normal;
						VectorD tangential_velocity = particles.u(i) - normal_velocity;
						VectorD displacement = normal * phi;
						particles.u(i) = tangential_velocity + displacement / dt;
					}
				}
			}
		}

		virtual void Advance(const int current_frame, const real current_time, const real dt) {
			Info("Advancing");
			Body_Force(dt);
			Update_Position(dt);
			if (Has_Boundary()) Enforce_Boundary(dt);
		}
	};
}