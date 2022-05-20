//////////////////////////////////////////////////////////////////////////
// Fluid Euler
// Copyright (c) (2022-), Bo Zhu, Mengdi Wang, Yitong Deng
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once

#include "Simulator.h"
#include "EulerParticles.h"
#include "MELPIO.h"
#include "SPH_Utils.h"

namespace Meso {
	template<int d>
	class MELP : public Simulator {
		Typedef_VectorD(d); Typedef_MatrixD(d);
	public:
		MELP() : Simulator() {}

		real enclosed_vol = 1.;
		real enclosed_amount = 1.;
		real p_out = 1.; 
		
		SPH<d> sph = SPH<d>(KernelType::QUINTIC);

		EulerParticles<d> e_particles;

		void Init() {
			e_particles.Update_Searcher();
			e_particles.Update_Local_Frames();
			e_particles.Orient_Normals();
			Update_Geometry();
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

		std::function<real(const int)> z_funcs(const int i) {
			std::function<real(const int)> z_func_i = [&](const int idx)->real {
				return (e_particles.x(idx) - e_particles.x(i)).dot(e_particles.Normal(i));
			};
			return z_func_i;
		}

		void Update_Metrics(void) {
			if constexpr (d == 3) {
#pragma omp parallel for
				for (int i = 0; i < e_particles.Size(); i++) {
					VectorT dzdx = sph.Gradient_Diff(e_particles.x(i), e_particles.E(i),
						z_funcs(i), 0., e_particles.aRef(), e_particles.xRef(), nbs, e_particles.Radius(i));
					MatrixT mt;
					real g11 = 1 + pow(dzdx[0], 2);
					real g22 = 1 + pow(dzdx[1], 2);
					real g12 = dzdx[0] * dzdx[1];
					mt << g11, g12, g12, g22;
					e_particles.g(i) = mt;
				}
			}
			else {
				e_particles.g(i) = MatrixT::Identity();
			}
		}

		VectorD Surface_Gradient(const int& i, std::function<real(const int)>& f) {
			Array<int> nbs; e_particles.nbs_searcher->Find_Nbs(e_particles.x(i), e_particles.Radius(i), nbs);
			VectorT tang_grad = sph.Gradient_Diff(e_particles.x(i), e_particles.E(i), f, f(i), e_particles.aRef(), e_particles.xRef(), nbs, e_particles.Radius(i));
			Vector2 tmp = e_particles.g(i) * tang_grad;
			Vector3 v; Unproject_To_World(tmp, e_particles.E(i), v);
			return v;
		}

		real Surface_Divergence(const int& i, std::function<VectorD(const int)>& f) {
			Array<int> nbs; e_particles.nbs_searcher->Find_Nbs(e_particles.x(i), e_particles.Radius(i), nbs);
			std::function<VectorT(const int)> project_orient_vec = [&](const int idx)->real {
				return Project_To_TPlane(f(idx)-f(i), e_particles.E(i));
			};
			MatrixT gg = e_particles.g(i).inverse();
			real sqrt_gg = sqrt(gg.determinant());
			real result = 0;
			for (int dim = 0; dim < d - 1; dim++) {
				std::function<real(const int)> tmp_func = [&](const int idx)->real {
					return sqrt_gg * project_orient_u(idx)[dim];
				};
				result += 1./sqrt_gg * sph.Gradient_Diff(e_particles.x(i), e_particles.E(i), tmp_func, tmp_func(i), e_particles.aRef(), e_particles.xRef(), nbs, e_particles.Radius(i));
			}
			return result;
		}

		void Update_Geometry(void) {
#pragma omp parallel for
			for (int i = 0; i < e_particles.Size(); i++) {
				Array<int> nbs;
				e_particles.nbs_searcher->Find_Nbs(e_particles.x(i), e_particles.Radius(i), nbs);
				e_particles.nden(i) = sph.Sum<real>(e_particles.x(i), e_particles.E(i),
					Array<real>(e_particles.Size(), 1.), e_particles.xRef(), nbs, e_particles.Radius(i));
				e_particles.a(i) = 1./e_particles.nden(i);
			}
#pragma omp parallel for
			for (int i = 0; i < e_particles.Size(); i++) {
				Array<int> nbs;
				e_particles.nbs_searcher->Find_Nbs(e_particles.x(i), e_particles.Radius(i), nbs);
				std::function<real(const int)> z_func_i = [&](const int idx)->real {
					return (e_particles.x(idx)-e_particles.x(i)).dot(e_particles.Normal(i));
				};
				e_particles.H(i) = sph.Laplacian<real>(e_particles.x(i), e_particles.E(i),
					z_func_i, 0., e_particles.aRef(), e_particles.xRef(), nbs, e_particles.Radius(i));
			}
			//Info("Total area: {}", e_particles.Size() * ArrayFunc::Mean<real>(e_particles.aRef()));
			//Info("curvature 0: {}", e_particles.H(0));
			enclosed_vol = Compute_Enclosed_Volume();
		}

		real Compute_Enclosed_Volume(void) {
			VectorD origin = VectorD::Zero();
			real vol = 0.;
			for (int i = 0; i < e_particles.Size(); i++) {
				real height = fabs(e_particles.Normal(i).dot(e_particles.x(i) - origin)); // volume of the skewed cone with the origin
				real cone_vol = real(1. / d) * e_particles.a(i) * height;
				if ((e_particles.x(i) - origin).dot(e_particles.Normal(i)) > 0.) {
					vol += cone_vol;
				}
				else {
					vol -= cone_vol;
				}
			}
			return vol;
		}

		void Update_Dynamics(const real dt) {
			real pressure = (enclosed_vol > 1.e-8) ? (enclosed_amount) / enclosed_vol : p_out;
			real pressure_diff = pressure - p_out;
			Info("Pressure_diff: {}", pressure_diff);
#pragma omp parallel for
			for (int i = 0; i < e_particles.Size(); i++) {
				e_particles.u(i) += 100 * dt * e_particles.Normal(i) * e_particles.H(i);
				e_particles.u(i) += 50 * dt * e_particles.Normal(i) * pressure_diff;
			}
			//tang velocity
			if constexpr (d == 3) {
#pragma omp parallel for
				for (int i = 0; i < e_particles.Size(); i++) {
					Array<int> nbs;
					e_particles.nbs_searcher->Find_Nbs(e_particles.x(i), e_particles.Radius(i), nbs);
					VectorT acc_2 = Surface_Gradient_Diff_SPH<d>(e_particles.x(i), e_particles.E(i),
						e_particles.ndenRef(), e_particles.nden(i), e_particles.aRef(), e_particles.xRef(), nbs,
						kernel, e_particles.Radius(i));
					VectorD acc_3; Unproject_To_World(acc_2, e_particles.E(i), acc_3);
					e_particles.u(i) += -30 * dt * acc_3;
				}
			}
		}

		virtual void Advance(const int current_frame, const real current_time, const real dt) {
			Info("Got here");
			Update_Geometry();
			Info("Got here 1");
			Update_Dynamics(dt);
#pragma omp parallel for
			for (int i = 0; i < e_particles.Size(); i++) {
				e_particles.x(i) += e_particles.u(i) * dt;
			}
		}
	};
}