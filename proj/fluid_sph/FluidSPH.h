//////////////////////////////////////////////////////////////////////////
// FluidSPH
// Copyright (c) (2022-), Yitong Deng
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once

#include "Simulator.h"
#include "SPHParticles.h"
#include "SPHIO.h"
#include "SPH_Utils.h"
#include "AnalyticalBoundary.h"

namespace Meso {

	template<int d>
	class FluidSPH : public Simulator {
		Typedef_VectorD(d); Typedef_MatrixD(d);
	public:
		FluidSPH() : Simulator() {}

		SPHParticles<d> particles;
		SPH_Utils<d> sph;
		AnalyticalBoundary<d> boundary;

		void Init() {
			Register_Nbs();
			//Test	
//			VectorD com = ArrayFunc::Mean<VectorD>(particles.xRef());
//#pragma omp parallel for
//			for (int i = 0; i < particles.Size(); i++) {
//				particles.u(i) = com - particles.x(i);
//				particles.u(i) *= 0.1;
//			}
		}

		virtual real CFL_Time(const real cfl) {
			real dx = particles.dx;
			real max_vel = ArrayFunc::Largest_Norm<VectorD>(particles.uRef());
			return dx * cfl / max_vel;
		}
		virtual void Output(DriverMetaData& metadata) {
			std::string vts_name = fmt::format("vtp{:04d}.vtp", metadata.current_frame);
			bf::path vtk_path = metadata.base_path / bf::path(vts_name);
			Write_SPH<d>(particles, vtk_path.string());
		}

		real Radius(const int i) const {
			return 4.0 * particles.dx;
		}

		int Boundary(const int i) const {
			return particles.B(i);
		}

		int Has_Boundary(void) const {
			return boundary.Available();
		}
			
		void Register_Nbs(void) {
			particles.Update_Searcher();
			if (particles.nbs_info.size() != particles.Size()) particles.nbs_info.resize(particles.Size());
#pragma omp parallel for
			for (int i = 0; i < particles.Size(); i++) {
				real r = Radius(i);
				Array<int> nbs; particles.nbs_searcher->Find_Nbs(particles.x(i), r, nbs);
				particles.nbs_info[i].set(nbs, r);
			}
		}

		void Pressure_Force(void) {
#pragma omp parallel for
			for (int i = 0; i < particles.Size(); i++) {
				VectorD grad_p = sph.Grad([&](int idx)->real {
					if (Boundary(idx)) return particles.rho(i);
					else return particles.rho(idx);
					},
					[&](int idx)->VectorD {return particles.x(i) - particles.x(idx); },
					[&](int idx)->real {return particles.V(idx); },
					particles.nbs_info[i](), particles.nbs_info[i].r);
				particles.acc(i) += -1. * grad_p;
			}
		}

		void Viscous_Force(void) {
//#pragma omp parallel for
//			for (int i = 0; i < particles.Size(); i++) {
//				VectorD lap_u = sph.Laplacian([&](int idx)->VectorD {return particles.u(idx) - particles.u(i); },
//					[&](int idx)->VectorD {return particles.x(i) - particles.x(idx); },
//					[&](int idx)->real {return particles.V(idx); },
//					particles.nbs_info[i](), particles.nbs_info[i].r);
//				particles.acc(i) += 0.025 * lap_u;
//			}
		}

		void Body_Force(void) {
#pragma omp parallel for
			for (int i = 0; i < particles.Size(); i++) {
				particles.acc(i) += -3.0 * VectorD::Unit(0);
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

		void Compute_Density(void) {
#pragma omp parallel for
			for (int i = 0; i < particles.Size(); i++) {
				particles.rho(i) = sph.Sum([&](int idx)->real {return particles.m(idx); },
					[&](int idx)->VectorD {return particles.x(i)-particles.x(idx); },
					particles.nbs_info[i](), particles.nbs_info[i].r);
				particles.V(i) = particles.m(i) / particles.rho(i);
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
						VectorD normal = particles.bnd_n(i);
						VectorD normal_velocity = particles.u(i).dot(normal) * normal;
						VectorD tangential_velocity = particles.u(i) - normal_velocity;
						VectorD displacement = normal * phi;
						particles.x(i) += displacement;
						particles.u(i) = tangential_velocity + 0.2 * displacement / dt;
					}
				}
			}
		}

		virtual void Advance(DriverMetaData& metadata) {
			Info("Advancing");
			real dt = metadata.dt;
			Prepare_Advance();
			Compute_Density();
			Pressure_Force();
			Viscous_Force();
			Body_Force();
			Update_Position(dt);
			if (Has_Boundary()) Enforce_Boundary(dt);
			Register_Nbs();
		}
	};
}