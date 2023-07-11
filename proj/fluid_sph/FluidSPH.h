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
		SPHParticles<d> b_particles;
		AnalyticalBoundary<d> boundary;
		Kernel kernel;

		void Init() {
			Register_Nbs(true);
		}

		virtual real CFL_Time(const real cfl) {
			real dx = particles.dx;
			real max_vel = ArrayFunc::Largest_Norm<VectorD, HOST>(particles.uRef());
			return dx * cfl / max_vel;
		}
		virtual void Output(DriverMetaData& metadata) {
			std::string vts_name = fmt::format("vtp{:04d}.vtp", metadata.current_frame);
			bf::path vtk_path = metadata.base_path / bf::path(vts_name);
			Write_SPH<d>(particles, b_particles, vtk_path.string());
		}

		real Radius(const int i) const {
			return 1.8 * particles.dx;
		}

		int Has_Boundary(void) const {
			return boundary.Available();
		}
			
		void Register_Nbs(bool init = false) {
			if (init) b_particles.Update_Searcher();
			particles.Update_Searcher();
			if (particles.nbs_info.size() != particles.Size()) particles.nbs_info.resize(particles.Size());
#pragma omp parallel for
			for (int i = 0; i < particles.Size(); i++) {
				real r = Radius(i);
				Array<int> nbs; particles.nbs_searcher->Find_Nbs(particles.x(i), r, nbs);
				Array<int> b_nbs; b_particles.nbs_searcher->Find_Nbs(particles.x(i), r, b_nbs);
				particles.nbs_info[i].set(nbs, b_nbs, r);
			}
		}

		void Compute_Pressure_WCSPH(void) {
			real rho_0 = 1.;
#pragma omp parallel for
			for (int i = 0; i < particles.Size(); i++) {
				particles.p(i) = 1. * (MathFunc::Power3(particles.rho(i) / rho_0) - 1.);
				particles.p(i) = std::max<real>(0., particles.p(i));
			}
		}

		void Compute_Pressure_IISPH(const real dt) {
			real rho_0 = 1.;
			real max_vel = ArrayFunc::Largest_Norm<VectorD>(particles.uRef());
			Info("max vel: {}", max_vel);
			Info("dt is: {}", dt);
			if (dt < 1.e-8) return; // if dt is 0. then return. Otherwise would be numerically troublesome
			real omega = 0.1; //Jacobi Relaxation
			int max_iter = 100;
			Array<real> a_ii(particles.Size(), 0.); //diagonal elements of matrix A
			Array<real> s_i(particles.Size(), 0.); //source term
			Array<VectorD> grad_p(particles.Size(), VectorD::Zero()); //gradient of p
			Array<real> LHS(particles.Size(), 0.); //Left hand side
			Array<real> err(particles.Size(), 0.); //error
			real avg_err = 0.;

			// init pressure to 0
#pragma omp parallel for
			for (int i = 0; i < particles.Size(); i++) {
				particles.p(i) = 0.;
			}

			// compute source
#pragma omp parallel for
			for (int i = 0; i < particles.Size(); i++) {
				real vel_div = sph.Div([&](int idx)->VectorD {return particles.u(idx) - particles.u(i); },
					[&](int idx)->VectorD {return particles.x(i) - particles.x(idx); },
					[&](int idx)->real {return particles.V(idx); },
					particles.nbs_info[i].nbs, particles.nbs_info[i].r);
				real pred_next_rho = particles.rho(i) + (dt * -1 * particles.rho(i) * vel_div);
				s_i[i] = rho_0 - pred_next_rho;
			}

#pragma omp parallel for
			for (int i = 0; i < particles.Size(); i++) {
				a_ii[i] = sph.Lap_Diagonal_Diff(particles.V(i),
					[&](int idx)->VectorD {return particles.x(i) - particles.x(idx); },
					[&](int idx)->real {return particles.V(idx); },
					particles.nbs_info[i].nbs, particles.nbs_info[i].r);
				//Info("a_ii: {}", a_ii[i]);
				a_ii[i] *= dt * dt;
				//Info("a_ii after: {}", a_ii[i]);
			}

			// Iterate
			for (int l = 0; l < max_iter; l++) {
#pragma omp parallel for
				for (int i = 0; i < particles.Size(); i++) { // for each E particle
					//Compute grad p
					grad_p[i] = sph.Grad([&](int idx)->real {return particles.p(idx) - particles.p(i); },
						[&](int idx)->VectorD {return particles.x(i) - particles.x(idx); },
						[&](int idx)->real {return particles.V(idx); },
						particles.nbs_info[i].nbs, particles.nbs_info[i].r);
				}
#pragma omp parallel for
				for (int i = 0; i < particles.Size(); i++) { // for each E particle
					//Compute LHS
					// second order
					real lap_p = sph.Div([&](int idx)->VectorD {return grad_p[idx] - grad_p[i]; },
						[&](int idx)->VectorD {return particles.x(i) - particles.x(idx); },
						[&](int idx)->real {return particles.V(idx); },
						particles.nbs_info[i].nbs, particles.nbs_info[i].r);
					LHS[i] = dt * dt * lap_p;
					//Update pressure
					particles.p(i) = std::max<real>(0., particles.p(i) + omega / a_ii[i] * (s_i[i] - LHS[i]));
					err[i] = abs(LHS[i] - s_i[i]);
				}
				real new_avg_err = ArrayFunc::Mean<real>(err);

				avg_err = new_avg_err;
				//Info("avg_err: {}", avg_err);
			}
		}

		void Compute_Pressure(const real dt) {
			Compute_Pressure_IISPH(dt);
			//Compute_Pressure_WCSPH();
		}

		void Pressure_Force(const real dt) {
			Compute_Pressure(dt);
#pragma omp parallel for
			for (int i = 0; i < particles.Size(); i++) {
				VectorD grad_p = sph.Grad([&](int idx)->real {return particles.p(idx) - particles.p(i); },
					[&](int idx)->VectorD {return particles.x(i) - particles.x(idx); }, 
					[&](int idx)->real {return particles.V(idx); }, 
					particles.nbs_info[i].nbs, particles.nbs_info[i].r);

				particles.acc(i) += -1./particles.rho(i) * grad_p;
			}
		}

		void Viscous_Force(void) {
#pragma omp parallel for
			for (int i = 0; i < particles.Size(); i++) {
				VectorD lap_u = sph.Lap_Brookshaw([&](int idx)->VectorD {return particles.u(idx) - particles.u(i); },
					[&](int idx)->VectorD {return particles.x(i) - particles.x(idx); },
					[&](int idx)->real {return particles.V(idx); },
					particles.nbs_info[i].nbs, particles.nbs_info[i].r);
				particles.acc(i) += 1. * lap_u;
			}
		}

		void Body_Force(void) {
#pragma omp parallel for
			for (int i = 0; i < particles.Size(); i++) {
				particles.acc(i) += -3.0 * VectorD::Unit(1);
			}
		}

		void Update_Position(const real dt) {
#pragma omp parallel for
			for (int i = 0; i < particles.Size(); i++) {
				particles.u(i) += dt * particles.acc(i);
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
					particles.nbs_info[i].nbs, particles.nbs_info[i].r);
				particles.rho(i) += sph.Sum([&](int idx)->real {return particles.m(i); },
					[&](int idx)->VectorD {return particles.x(i) - b_particles.x(idx); },
					particles.nbs_info[i].b_nbs, particles.nbs_info[i].r);
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

		virtual void Advance(DriverMetaData& metadata) {
			Info("Advancing");
			real dt = metadata.dt;
			Prepare_Advance();
			Compute_Density();
			Pressure_Force(dt);
			Viscous_Force();
			Body_Force();
			Update_Position(dt);
			if (Has_Boundary()) Enforce_Boundary(dt);
			Register_Nbs();
		}
	};
}