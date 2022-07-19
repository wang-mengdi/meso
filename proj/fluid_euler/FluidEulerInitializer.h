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
	};
}