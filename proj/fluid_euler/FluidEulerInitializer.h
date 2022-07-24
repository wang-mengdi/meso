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
#include "FluidClebsch.h"
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

			/*Grid<d> sol_grid = sol.grid;
			sol.Iterate_Nodes(
				[&](const VectorDi cell) {
					std::cout << "cell is " << cell[0] <<" " << cell[1];
					for (int i = 0; i < 4; i++) {
						VectorDi nei_cell = sol_grid.Neighbor_Node(cell, i);
						if (sol_grid.Valid(nei_cell)) {
							int nei_index = sol_grid.Index(nei_cell);
							std::cout << " neighbor " << i << ", fixed is " << fixed.Data()[nei_index] << ", vol is " << vol.Data(0)[nei_index] << " " << vol.Data(1)[nei_index];
						}
					}
					std::cout << std::endl;
				}
			);*/
			Field<real, d> sol_host = sol;

			// read sol_complex from complex_sol.txt
			/*FieldDv<real, d> sol_complex(grid);
			
			std::ifstream in("complex_sol.txt");
			std::string line;
			std::getline(in, line);
			in.close();

			std::stringstream inStream(line);
			std::string each;

			for (int i = 0; i != 16; i++) {
				for (int j = 0; j != 16; j++) {
					real temp;
					VectorDi cell;
					cell[0] = i;
					cell[1] = j;
					std::getline(inStream, each, ',');
					std::stringstream tmpStream(each);
					tmpStream >> temp;
					sol_complex.Data()[sol_complex.grid.Index(cell)] = temp;
				}
			}

			Info("complex_sol is {}", sol_complex);
			//Info("Start iterate nodes");

			FieldDv<real, d> sol_complex_copy = sol_complex;
			//auto [iter_comp, error_comp] = MGPCG.Solve(sol_complex_copy.Data(), wz.Data());
			//Info("iter_comp is {}, error_comp is {}", iter_comp, error_comp);
			poisson.Apply(sol_complex_copy.Data(), sol_complex.Data());
			sol_complex_copy -= wz;
			Info("Residual is {}", sol_complex_copy);

			Info("sol is {}", sol);
			Info("grid.dx is {}", grid.dx);*/
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
			//Info("initial_vel is {}", initial_vel);
		}

		// trefoil knot
		void Case_4(json& j, FluidEuler<d>& fluid) {
			FluidClebsch<d> fluid_clebsch;
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
			Initialize_Trefoil_Knot_Complex_Poly(grid, initial_complex_poly, center, radius);
			Info("initial_complex_poly is {}", initial_complex_poly);

			real h_bar = Json::Value(j, "h_bar", 0.01);

			// transform initial complex poly value into psi
			//Transform_Complex_Poly_To_Psi(h_bar, grid, initial_complex_poly, initial_wave_func);
			initial_wave_func = initial_complex_poly;

			fluid_clebsch.Init(h_bar, fixed, initial_wave_func, vol);
			//Vector2C a = { C(1, 0), C(0.5, 0.5) };
			//Vector2C b = { C(1, 0), C(0, 2.0) };
			//ArrayFunc::Conj_Dot(a, b);
			initial_vel = fluid_clebsch.velocity;
			fluid.Init(fixed, vol, face_fixed, initial_vel);
			Info("initial_vel is {}", initial_vel);
		}

		void Initialize_Trefoil_Knot_Complex_Poly(const Grid<d> grid, Field<Vector2C, d>& initial_complex_poly, const VectorD& center, const VectorD& radius) {
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

					initial_complex_poly(cell) = ArrayFunc::V2C(ArrayFunc::C2V(psi).normalized());
				}
			);
		}

		//void Transform_Complex_Poly_To_Psi(real h_bar, const Grid<d> grid, Field<Vector4, d>& initial_complex_poly, Field<Vector2C, d>& initial_wave_func) {
		//	initial_wave_func.Init(grid);

		//	Eigen::Matrix<int, 3, 2> bc_width;
		//	Eigen::Matrix<real, 3, 2> bc_val;
		//	bc_width << -1, -1, -1, -1, -1, -1;
		//	bc_val << 0, 0, 0, 0, 0, 0;

		//	std::shared_ptr<Array<Vector4, DEVICE>> data1 = std::make_shared<Array<Vector4, DEVICE>>(grid.Memory_Size());
		//	std::shared_ptr<Array<Vector4, DEVICE>> data2 = std::make_shared<Array<Vector4, DEVICE>>(grid.Memory_Size());
		//	ArrayFunc::Fill(*data1, Vector4(0));
		//	ArrayFunc::Fill(*data2, Vector4(0));
		//	ArrayFunc::Dot(*data1, *data2);
		//	MaskedPoissonMapping<Vector4, d> poisson;
			//VCycleMultigridIntp<Vector4, d> MG_precond;
			//ConjugateGradient<Vector4> MGPCG;
			//Field<Vector4, d> P_Q_Laplace(grid);
			//FaceField<real, d> alpha(grid.Counts(), (real)1);
			//// question: Can I directly use this function when data in field is Vector2C
			//// if not : if I change Vector2C to Vector4, will it work?
			//poisson.Init(fixed, alpha);
			//poisson.Apply(P_Q_Laplace.Data(), initial_complex_poly.Data());

			//// compute rhs of poisson equation which is used to compute q
			//Field<real, d> q_rhs_host(grid);
			//grid.Exec_Nodes(
			//	[&](const Vector<int, d> cell) {
			//		Vector<C, 2> psi_tmp_1 = V2C(P_Q_Laplace(cell));
			//		Vector<C, 2> psi_tmp_2 = V2C(initial_complex_poly(cell));
			//		C q_rhs_tmp = ArrayFunc::Conj_Dot(psi_tmp_1, psi_tmp_2);
			//		q_rhs_host(cell) = q_rhs_tmp.real() * h_bar;
			//	}
			//);

			//// solve q
			//FieldDv<real, d> q_rhs;
			//q_rhs = q_rhs_host;
			//// question: Do I need this?  q_rhs *= dx * dx; 

			//FieldDv<real, d> q(grid.Counts(), (real)0);

			//poisson.Init(fixed, alpha);
			//MG_precond.Init_Poisson(poisson, 2, 2);
			//MGPCG.Init(&poisson, &MG_precond, false, -1, 1e-6);

			//auto [iter, error] = MGPCG.Solve(q.Data(), q_rhs.Data());
			//Info("iter is {}, error is {}", iter, error);

			//// using q and complex poly value to calculate psi
			//grid.Exec_Nodes(
			//	[&](const Vector<int, d> cell) {
			//		const VectorD pos = grid.Position(cell);
			//		Vector2C psi = initial_complex_poly(cell);
			//		for (int i = 0; i < 2; i++) psi[i] *= C(thrust::exp(-C(0., 1.) * q(cell) / h_bar));
			//		initial_wave_func(cell) = psi;
			//	}
			//);
		//}




	};
}