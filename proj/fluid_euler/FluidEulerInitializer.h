//////////////////////////////////////////////////////////////////////////
// Initializer of a Fluid Euler System
// Copyright (c) (2022-), Bo Zhu, Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once

#include "FluidEuler.h"
#include "ImplicitManifold.h"
#include "Json.h"
#include "GridEulerFunc.h"
#include "WaveFunc.h"
#include <fstream>
#include <iostream>
#include <sstream>


namespace Meso {
	template<int d>
	class  FluidEulerInitializer{
	public:
		Typedef_VectorD(d);

		//define the boundary conditions of the eulerian system
		Field<bool, d> fixed;
		FaceField<real, d> vol;
		FaceField<bool, d> face_fixed;
		FaceField<real, d> initial_vel;

		void Apply(json& j, FluidEuler<d>& fluid) {
			int test = Json::Value(j, "test", 0);
			switch (test) {
			case 0:Case_0(j, fluid); break;
			case 1:Case_1(j, fluid); break;
			case 2:Case_2(j, fluid); break;
			case 3:Case_3(j, fluid); break;
			case 4:Case_4(j, fluid); break;
			case 5:Case_5(j, fluid); break;
			default:Assert(false, "test {} not exist", test); break;
			}
		}

		//an empty field with inflow boundary condition v=1
		void Case_0(json& j, FluidEuler<d>& fluid) {
			int scale = Json::Value(j, "scale", 32);
			real dx = 1.0 / scale;
			VectorDi grid_size = scale * MathFunc::Vi<d>(2, 1, 1);
			Grid<d> grid(grid_size, dx, VectorD::Zero(), CENTER);

			Eigen::Matrix<int, 3, 2> bc_width;
			Eigen::Matrix<real, 3, 2> bc_val;
			bc_width << 0, -1, 0, 0, 0, 0;
			bc_val << 1, 1, 0, 0, 0, 0;

			GridEulerFunc::Set_Boundary(grid, bc_width, bc_val, fixed, vol, face_fixed, initial_vel);
			fluid.Init(fixed, vol, face_fixed, initial_vel);
			ArrayFunc::Fill(fluid.velocity.Data(0), 1.0);
		}

		void Case_1(json& j, FluidEuler<d>& fluid) {
			int scale = Json::Value(j, "scale", 32);
			real dx = 1.0 / scale;
			VectorDi grid_size = scale * MathFunc::Vi<d>(1, 1, 1);
			Grid<d> grid(grid_size, dx, VectorD::Zero(), CENTER);

			Eigen::Matrix<int, 3, 2> bc_width;
			Eigen::Matrix<real, 3, 2> bc_val;
			bc_width << 1, 1, 1, 1, 1, 1;
			bc_val << 0, 0, 0, 0, 0, 0;

			GridEulerFunc::Set_Boundary(grid, bc_width, bc_val, fixed, vol, face_fixed, initial_vel);

			grid.Exec_Faces(
				[&](const int axis, const VectorDi face) {
					if (axis == 0 && face[1] >= grid.Counts()[1] - 2) {
						if (!face_fixed(axis, face)) {
							face_fixed(axis, face) = true;
							initial_vel(axis, face) = 1.0;
						}
					}
				}
			);

			fluid.Init(fixed, vol, face_fixed, initial_vel);
		}

		// with a sheprical obstacle, initial velocity of 1
		void Case_2(json& j, FluidEuler<d>& fluid) {
			int scale = Json::Value(j, "scale", 32);
			real dx = 1.0 / scale;
			VectorDi grid_size = scale * MathFunc::Vi<d>(2, 1, 1);
			Grid<d> grid(grid_size, dx, VectorD::Zero(), CENTER);

			Eigen::Matrix<int, 3, 2> bc_width;
			Eigen::Matrix<real, 3, 2> bc_val;
			bc_width << 0, 0, 0, 0, 0, 0;
			bc_val << 1, 1, 0, 0, 0, 0;

			////sphere
			VectorD center = grid.Center(); 
			center[0] = center[0] / 2;
			real r = (real)0.1;
			Sphere<d> sphere(center, r);

			GridEulerFunc::Set_Boundary(grid, bc_width, bc_val, fixed, vol, face_fixed, initial_vel);

			grid.Exec_Nodes(
				[&](const VectorDi cell) {
					const VectorD pos = grid.Position(cell);
					if (sphere.Inside(pos)) {
						fixed(cell) = true;
					}
				}
			);

			grid.Exec_Faces(
				[&](const int axis, const VectorDi face) {
					const VectorD pos = grid.Face_Center(axis, face);
					if (sphere.Inside(pos)) {
						face_fixed(axis, face) = true;
						initial_vel(axis, face) = 0.0;
						vol(axis, face) = 0;
					}
				}
			);
			fluid.Init(fixed, vol, face_fixed, initial_vel);
		}


		// taylor vortex
		void Case_3(json& j, FluidEuler<d>& fluid) {
			int scale = Json::Value(j, "scale", 32);
			const real pi = 3.1415926;
			real dx = 2.0 * pi / scale;
			VectorDi grid_size = scale * MathFunc::Vi<d>(1, 1, 1);
			Grid<d> grid(grid_size, dx, VectorD::Zero(), CENTER);
			MaskedPoissonMapping<real, d> poisson;
			VCycleMultigridIntp<real, d> MG_precond;
			ConjugateGradient<real> MGPCG;

			Eigen::Matrix<int, 3, 2> bc_width;
			Eigen::Matrix<real, 3, 2> bc_val;
			bc_width << -1, -1, -1, -1, -1, -1;
			bc_val << 0, 0, 0, 0, 0, 0;

			GridEulerFunc::Set_Boundary(grid, bc_width, bc_val, fixed, vol, face_fixed, initial_vel);
			

			VectorD vortex_p1, vortex_p2;
			// two vortices are 0.81 apart
			vortex_p1[0] = grid.Center()[0] + (real)0.405;
			vortex_p2[0] = grid.Center()[0] - (real)0.405;
			vortex_p1[1] = grid.Center()[1]; 
			vortex_p2[1] = vortex_p1[1];

			Field<real, d> wz_host(grid);

			grid.Exec_Nodes(
				[&](const VectorDi cell) {
					const VectorD pos = grid.Position(cell);
					real pr1 = pow((pos[0] - vortex_p1[0]), 2) + pow((pos[1] - vortex_p1[1]), 2);
					real pr2 = pow((pos[0] - vortex_p2[0]), 2) + pow((pos[1] - vortex_p2[1]), 2);
					wz_host(cell) = (real)1 / (real)0.3 * ((real)2 - pr1 / (real)0.09) * exp((real)0.5 * ((real)1 - pr1 / (real)0.09));
					wz_host(cell) += (real)1 / (real)0.3 * ((real)2 - pr2 / (real)0.09) * exp((real)0.5 * ((real)1 - pr2 / (real)0.09));
					wz_host(cell) *= -grid.dx * grid.dx;
				}
			);
			FieldDv<real, d> wz;
			wz = wz_host;
			//Info("wz is {}", wz);

			FieldDv<real, d> sol(grid.Counts(), (real)0);
			FaceField<real, d> alpha(grid.Counts(), (real)1);

			poisson.Init(fixed, alpha);
			MG_precond.Init_Poisson(poisson, 2, 2);
			MGPCG.Init(&poisson, &MG_precond, false, -1, 1e-6);

			auto [iter, error] = MGPCG.Solve(sol.Data(), wz.Data());
			Info("iter is {}, error is {}", iter, error);

			Field<real, d> sol_host = sol;

			grid.Exec_Faces(
				[&](const int axis, const VectorDi face) {
					const VectorD pos = grid.Face_Center(axis, face);

					if (face[axis] == 0 || face[axis] == grid.Counts()[axis]) { return; }

					const VectorDi cell_left = face - VectorDi::Unit(axis);
					const VectorDi cell_right = face;

					if (axis == 0) {
						initial_vel(1, face) = (sol_host(cell_right) - sol_host(cell_left)) / grid.dx;
					}
					else if (axis == 1) {
						initial_vel(0, face) = (sol_host(cell_left) - sol_host(cell_right)) / grid.dx;
					}
				}
			);

			fluid.Init(fixed, vol, face_fixed, initial_vel);
		}

		// trefoil knot
		void Case_4(json& j, FluidEuler<d>& fluid) {
			Field<Vector2C, d> initial_wave_func;
			int scale = Json::Value(j, "scale", 32);
			const real pi = acos(-1.0);
			real dx = (real)2.0 * pi / scale;
			VectorDi grid_size = scale * MathFunc::Vi<d>(1, 1, 1);
			Grid<d> grid(grid_size, dx, VectorD::Zero(), CENTER);

			Eigen::Matrix<int, 3, 2> bc_width;
			Eigen::Matrix<real, 3, 2> bc_val;
			bc_width << -1, -1, -1, -1, -1, -1;
			bc_val << 0, 0, 0, 0, 0, 0;

			GridEulerFunc::Set_Boundary(grid, bc_width, bc_val, fixed, vol, face_fixed, initial_vel);

			Field<Vector2C, d> initial_complex_poly;
			VectorD center = MathFunc::V<d>(pi, pi, 0.6 * pi);
			VectorD radius = MathFunc::V<d>(1.2, 1.2, 1.6);
			initial_complex_poly.Init(grid);
			grid.Exec_Nodes(
				[&](const Vector<int, d> cell) {
					const VectorD pos = grid.Position(cell);
					VectorD XYZ = pos - center; for (int i = 0; i < d; i++) { XYZ[i] *= radius[i]; }
					real R = XYZ.norm();

					real fR = exp(-pow(R, 8) / pow(9, 4));
					C alpha = C(XYZ[0], XYZ[1]) * 2. * fR / (1 + R * R);
					//C beta = C(2.*(XYZ[2]-1)*fR, (1+R*R))/(1+R*R);
					C beta = C(2. * XYZ[2] * fR, (1 + R * R - 2 * fR)) / (1 + R * R);
					C P = alpha * alpha * alpha * alpha * alpha;
					C Q = alpha * alpha * alpha + beta * beta;
					//P: psi1; Q: psi2
					Vector<C, 2> psi;
					psi[0] = P;
					psi[1] = Q;

					initial_complex_poly(cell) = WaveFunc::Vector4_To_Vector2C(WaveFunc::Vector2C_To_Vector4(psi).normalized());
				}
			);

			real h_bar = Json::Value(j, "h_bar", 0.01);

			// transform initial complex poly value into psi
			initial_wave_func = initial_complex_poly;

			FaceFieldDv<real, d> initial_vel_dev = initial_vel;
			FieldDv<Vector2C, d> initial_wave_func_dev = initial_wave_func;
			WaveFunc::Exterior_Derivative_W2V(initial_vel_dev, initial_wave_func_dev, h_bar);
			initial_vel = initial_vel_dev;
			fluid.Init(fixed, vol, face_fixed, initial_vel);
		}

		// leapfrog
		void Case_5(json& j, FluidEuler<d>& fluid) {
			Field<Vector2C, d> initial_wave_func;
			int scale = Json::Value(j, "scale", 32);

			real pi = 3.1415926;
			real length = 2 * pi;
			real dx = length / scale;
			VectorDi grid_size = scale * MathFunc::Vi<d>(1, 1, 1);
			Grid<d> grid(grid_size, dx, VectorD::Zero(), CENTER);

			Eigen::Matrix<int, 3, 2> bc_width;
			Eigen::Matrix<real, 3, 2> bc_val;
			bc_width << -1, -1, -1, -1, -1, -1;
			bc_val << 0, 0, 0, 0, 0, 0;

			GridEulerFunc::Set_Boundary(grid, bc_width, bc_val, fixed, vol, face_fixed, initial_vel);

			initial_wave_func.Init(grid);
			real background_speed=(real)-.0;//2;
			VectorD bdry_v=VectorD::Unit(0)*(real)-.0;//2;

			////initialize leapfrogging velocity
			const real h_bar = 0.5;
			grid.Exec_Nodes(
				[&](const Vector<int, d> cell) {
					const VectorD pos = grid.Position(cell);
					Vector2C psi = WaveFunc::Vel_To_Psi_C<d>(bdry_v, pos);
					Vector<VectorD, 2> c = { MathFunc::V<d>((real)length * .5,(real)length / 2,(real)length / 2), MathFunc::V<d>((real)length * .5,(real)length / 2,(real)length / 2) };
					Vector<real, 2> r = { 1.8,1. };
					int n = (int)c.size();
					for (int i = 0; i < n; i++) {
						real rx = (pos[0] - c[i][0]) / r[i];
						real r2 = (pos - c[i]).squaredNorm() / pow(r[i], 2);
						real D = exp(-pow(r2 / (real)9, (real)4));
						C q(2. * rx * D / (r2 + 1), (r2 + 1. - 2. * D) / (r2 + 1));
						psi[0] *= q;
					}
					initial_wave_func(cell) = WaveFunc::Vector4_To_Vector2C(WaveFunc::Vector2C_To_Vector4(psi).normalized());
				}
			);
			FaceFieldDv<real, d> initial_vel_dev = initial_vel;
			FieldDv<Vector2C, d> initial_wave_func_dev = initial_wave_func;
			WaveFunc::Exterior_Derivative_W2V(initial_vel_dev, initial_wave_func_dev, h_bar);
			initial_vel = initial_vel_dev;
			fluid.Init(fixed, vol, face_fixed, initial_vel);
		}

	};
}