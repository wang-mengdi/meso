//////////////////////////////////////////////////////////////////////////
// MELP
// Copyright (c) (2022-), Bo Zhu, Mengdi Wang, Yitong Deng
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once

#include "Simulator.h"
#include "EulerParticles.h"
#include "MELPIO.h"
#include "SPH.h"
#include "ParticleFrames.h"

namespace Meso {

	template<int d>
	class MELP : public Simulator {
		Typedef_VectorD(d); Typedef_MatrixD(d);
	public:
		MELP() : Simulator() {}

		real enclosed_vol = 1.;
		real enclosed_amount = 1.;
		real p_out = 1.; 
		
		SPH<d> sph;

		EulerParticles<d> e_particles;

		void Init() {
			Info("1");
			e_particles.Update_Searcher();
			Info("1.5");
			Register_Nbs();
			Info("2");
			e_particles.Update_Local_Frames();
			Info("3");
			e_particles.Orient_Normals();
			Info("4");
			Update_Geometry();
			Info("5");
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

		void Update_Control_Area(void) {
#pragma omp parallel for
			for (int i = 0; i < e_particles.Size(); i++) {
				const Array<int>& nbs = e_particles.nbs_info[i].e();
				e_particles.nden(i) = sph.template Sum<real>(e_particles.x(i), e_particles.E(i),
					Array<real>(e_particles.Size(), 1.), e_particles.xRef(), nbs, e_particles.Radius(i));
				e_particles.a(i) = 1. / e_particles.nden(i);
			}
			Info("Total area: {}", e_particles.Size() * ArrayFunc::Mean<real>(e_particles.aRef()));
			//Info("curvature 0: {}", e_particles.H(0));
		}

		void Update_Metric(void) {
			if constexpr (d == 3) {
#pragma omp parallel for
				for (int i = 0; i < e_particles.Size(); i++) {
					const Array<int>& nbs = e_particles.nbs_info[i].e();
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
#pragma omp parallel for
				for (int i = 0; i < e_particles.Size(); i++) {
					e_particles.g(i) = MatrixT::Identity();
				}
			}
		}

		void Update_Curvature(void) {
			if constexpr(d==3){
#pragma omp parallel for
				for (int i = 0; i < e_particles.Size(); i++) {
					const Array<int>& nbs = e_particles.nbs_info[i].e();
					std::function<real(const int)> z_func_i = [&](const int idx)->real {
						return (e_particles.x(idx) - e_particles.x(i)).dot(e_particles.Normal(i));
					};
					e_particles.H(i) = sph.template Laplacian<real>(e_particles.x(i), e_particles.E(i),
						z_func_i, 0., e_particles.aRef(), e_particles.xRef(), nbs, e_particles.Radius(i));
				}
				Info("Position 0 \n {}", e_particles.x(0));
				VectorD tmp_COM = VectorD::Zero();
				VectorD true_normal = (e_particles.x(0)- tmp_COM).normalized();
				VectorD true_curvature_vec = true_normal * 2. / 1.;
				Info("true curvature vec \n {}", -true_curvature_vec);
				Info("Curvature vector 0 \n {}", e_particles.H(0) * e_particles.Normal(0));
				Array<VectorD> curvature_normals(e_particles.Size());
				for (int dim = 0; dim < 3; dim++) {
					Array<VectorD> grad(e_particles.Size());
#pragma omp parallel for
					for (int i = 0; i < e_particles.Size(); i++) {
						const Array<int>& nbs = e_particles.nbs_info[i].e();
						std::function<real(const int)> x_func = [&](const int idx)->real {
							return e_particles.x(idx)[dim];
						};
						grad[i] = Surface_Gradient(i, x_func);
					}
#pragma omp parallel for
					for (int i = 0; i < e_particles.Size(); i++) {
						const Array<int>& nbs = e_particles.nbs_info[i].e();
						std::function<VectorD(const int)> grad_func = [&](const int idx)->VectorD {
							return grad[idx];
						};
						curvature_normals[i][dim] = Surface_Divergence(i, grad_func);
					}
				}
				Info("Laplace Beltrami vector 0 \n {}", curvature_normals[0]);
			}
		}

		VectorD Surface_Gradient(const int i, std::function<real(const int)>& f) {
			const Array<int>& nbs = e_particles.nbs_info[i].e();
			VectorT tang_grad = sph.Gradient_Diff(e_particles.x(i), e_particles.E(i), f, f(i), e_particles.aRef(), e_particles.xRef(), nbs, e_particles.Radius(i));
			Vector2 tmp = e_particles.g(i) * tang_grad;
			Vector3 v; Unproject_To_World(tmp, e_particles.E(i), v);
			return v;
		}

		VectorT TPlane(VectorD& u, MatrixD& E) {
			VectorT t_coords;
			for (int i = 0; i < d - 1; i++) {
				t_coords[i] = u.dot(E.col(i));
			}
			return t_coords;
		}

		real Surface_Divergence(const int i, std::function<VectorD(const int)>& f) {
			const Array<int>& nbs = e_particles.nbs_info[i].e();
			std::function<VectorT(const int)> project_orient_vec = [&](const int idx)->VectorT {
				VectorD thing = f(idx);
				return TPlane(thing, e_particles.E(i));
				//return VectorT::Zero();
			};
			MatrixT gg = e_particles.g(i).inverse();
			real sqrt_gg = sqrt(gg.determinant());
			real result = 0;
			for (int dim = 0; dim < d - 1; dim++) {
				std::function<real(const int)> tmp_func = [&](const int idx)->real {
					return sqrt_gg * project_orient_vec(idx)[dim];
				};
				real shit = sph.Gradient_Diff(e_particles.x(i), e_particles.E(i), tmp_func, tmp_func(i), e_particles.aRef(), e_particles.xRef(), nbs, e_particles.Radius(i))[dim];
				result += 1. / sqrt_gg * shit;
			}
			return result;
		}
		

		void Update_Geometry(void) {
			Info("A");
			Update_Control_Area();
			Info("B");
			Update_Metric();
			Info("C");
			Update_Curvature();
			Info("D");
//#pragma omp parallel for
//			for (int i = 0; i < e_particles.Size(); i++) {
//				Array<int> nbs;
//				e_particles.nbs_searcher->Find_Nbs(e_particles.x(i), e_particles.Radius(i), nbs);
//				e_particles.nden(i) = sph.Sum<real>(e_particles.x(i), e_particles.E(i),
//					Array<real>(e_particles.Size(), 1.), e_particles.xRef(), nbs, e_particles.Radius(i));
//				e_particles.a(i) = 1./e_particles.nden(i);
//			}
//#pragma omp parallel for
//			for (int i = 0; i < e_particles.Size(); i++) {
//				Array<int> nbs;
//				e_particles.nbs_searcher->Find_Nbs(e_particles.x(i), e_particles.Radius(i), nbs);
//				std::function<real(const int)> z_func_i = [&](const int idx)->real {
//					return (e_particles.x(idx)-e_particles.x(i)).dot(e_particles.Normal(i));
//				};
//				e_particles.H(i) = sph.Laplacian<real>(e_particles.x(i), e_particles.E(i),
//					z_func_i, 0., e_particles.aRef(), e_particles.xRef(), nbs, e_particles.Radius(i));
//			}
//			//Info("Total area: {}", e_particles.Size() * ArrayFunc::Mean<real>(e_particles.aRef()));
//			//Info("curvature 0: {}", e_particles.H(0));
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
			//return 0.;
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
					const Array<int>& nbs = e_particles.nbs_info[i].e();
					std::function<real(const int)> nden_func = [&](const int idx)->real {
						return e_particles.nden(idx);
					};
					VectorD acc_3 = Surface_Gradient(i, nden_func);
					e_particles.u(i) += -30 * dt * acc_3;
				}
			}
		}

		void Register_E_Nbs(void) {
			if (e_particles.nbs_info.size() != e_particles.Size()) e_particles.nbs_info.resize(e_particles.Size());
#pragma omp parallel for
			for (int i = 0; i < e_particles.Size(); i++) {
				Array<int> nbs; e_particles.nbs_searcher->Find_Nbs(e_particles.x(i), e_particles.Radius(i), nbs);
				e_particles.nbs_info[i].set_e(nbs);
			}
		}

		void Register_Nbs(void) {
			Register_E_Nbs();
		}

		virtual void Advance(const int current_frame, const real current_time, const real dt) {
			Register_Nbs();
			Info("num e nbs 0: {}", e_particles.nbs_info[0].size_e());
			Info("num l nbs 0: {}", e_particles.nbs_info[0].size_l());
			Info("Got here");
			//Update_Control_Area();
			//Update_Metrics();
			//Update_Curvature();
			std::cout << "metrics 0: \n" << e_particles.g(0) << std::endl;
			Update_Geometry();
//			Info("Got here 1");
			Update_Dynamics(dt);
#pragma omp parallel for
			for (int i = 0; i < e_particles.Size(); i++) {
				e_particles.x(i) += e_particles.u(i) * dt;
			}
		}
	};
}