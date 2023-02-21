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
#include "Random.h"
#include <fstream>
#include <iostream>
#include <sstream>


namespace Meso {
	template<int d>
	class  FluidEulerInitializer{
	public:
		Typedef_VectorD(d);

		//define the boundary conditions of the eulerian system
		Field<unsigned char, d> cell_type;
		FaceField<real, d> vol;
		FaceField<bool, d> face_fixed;
		FaceField<real, d> initial_vel;

		void Apply(json& j, FluidEuler<d>& fluid) {
			std::string test = Json::Value(j, "test", std::string("taylor"));
			if (test == "taylor")Init_Taylor(j, fluid);
			else if (test == "leapfrog")Init_Leapfrog(j, fluid);
			//else if (test == "karman")Init_Karman(j, fluid);
			else Assert(false, "test {} not exist", test);
		}

		// taylor vortex
		void Init_Taylor(json& j, FluidEuler<d>& fluid) {
			int scale = Json::Value(j, "scale", 32);
			const real pi = 3.1415926;
			real dx = 2.0 * pi / scale;
			VectorDi grid_size = scale * MathFunc::Vi<d>(1, 1, 1);
			Grid<d> grid(grid_size, dx, VectorD::Zero(), CENTER);
			MaskedPoissonMapping<real, d> poisson;
			VCycleMultigridIntp<real, d> MG_precond;
			ConjugateGradient<real, d> MGPCG;

			Eigen::Matrix<int, 3, 2> bc_width;
			Eigen::Matrix<real, 3, 2> bc_val;
			bc_width << -1, -1, -1, -1, -1, -1;
			bc_val << 0, 0, 0, 0, 0, 0;

			GridEulerFunc::Set_Boundary(grid, bc_width, bc_val, cell_type, vol, face_fixed, initial_vel);
			

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
			
			poisson.Init(cell_type, alpha);
			MG_precond.Init_Poisson(poisson, 2);
			MGPCG.Init(&poisson, &MG_precond, false, -1, 1e-5);

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

			fluid.Init(cell_type, vol, face_fixed, initial_vel);
		}


		//covector leapfrog
		void Init_Leapfrog(json& j, FluidEuler<d>& fluid) {
			int scale = Json::Value(j, "scale", 32);

			real pi = 3.1415926;
			real length = 2 * pi;
			real dx = length / scale;
			VectorDi grid_size = scale * MathFunc::Vi<d>(1, 1, 1);
			Grid<d> grid(grid_size, dx, VectorD::Zero(), CENTER);

			Eigen::Matrix<int, 3, 2> bc_width;
			Eigen::Matrix<real, 3, 2> bc_val;
			bc_width << 0, 0, 0, 0, 0, 0;
			bc_val << 0, 0, 0, 0, 0, 0;
			GridEulerFunc::Set_Boundary(grid, bc_width, bc_val, cell_type, vol, face_fixed, initial_vel);

			if constexpr (d == 2) {
				Vector2d center = grid.Center();
				Field<real, d> wz_host(grid);
				real a = 0.02;
				grid.Exec_Nodes(
					[&](const VectorDi cell) {
						Vector<real, 2> pos = VectorD(cell(0), cell(1)) * grid.dx - VectorD::Ones() * pi;
				Vector<real, 2> vort_pos0 = Vector<real, 2>(-0.75, -2);
				Vector<real, 2> vort_pos1 = Vector<real, 2>(0.75, -2);
				Vector<real, 2> vort_pos2 = Vector<real, 2>(-1.5, -2);
				Vector<real, 2> vort_pos3 = Vector<real, 2>(1.5, -2);
				real r_sqr0 = pow(pos[0] - vort_pos0[0], 2) + pow(pos[1] - vort_pos0[1], 2);
				real r_sqr1 = pow(pos[0] - vort_pos1[0], 2) + pow(pos[1] - vort_pos1[1], 2);
				real r_sqr2 = pow(pos[0] - vort_pos2[0], 2) + pow(pos[1] - vort_pos2[1], 2);
				real r_sqr3 = pow(pos[0] - vort_pos3[0], 2) + pow(pos[1] - vort_pos3[1], 2);
				real c0 = 1000 / (2 * pi) * exp(-0.5 * r_sqr0 / a / a);
				real c1 = -1000 / (2 * pi) * exp(-0.5 * r_sqr1 / a / a);
				real c2 = 1000 / (2 * pi) * exp(-0.5 * r_sqr2 / a / a);
				real c3 = -1000 / (2 * pi) * exp(-0.5 * r_sqr3 / a / a);
				wz_host(cell) = c0 + c1 + c2 + c3;
				wz_host(cell) *= grid.dx * grid.dx;
					});
				FieldDv<real, d> wz;
				wz = wz_host;

				FieldDv<real, d> sol(grid.Counts(), (real)0);
				FaceField<real, d> alpha(grid.Counts(), (real)1);

				MaskedPoissonMapping<real, d> poisson;
				VCycleMultigridIntp<real, d> MG_precond;
				ConjugateGradient<real, d> MGPCG;
				poisson.Init(cell_type, alpha);
				MG_precond.Init_Poisson(poisson, 2);
				MGPCG.Init(&poisson, &MG_precond, false, -1, 1e-6, true);

				auto [iter, error] = MGPCG.Solve(sol.Data(), wz.Data());
				Info("iter is {}, error is {}", iter, error);

				Field<real, d> sol_host = sol;
				Field<real, d> sol_host_big;
				Grid<d> big_grid(grid.Counts() + VectorDi::Ones());
				sol_host_big.Init(big_grid, 0);
				grid.Iterate_Nodes([&](const VectorDi cell)
					{
						sol_host_big(cell) = sol_host(cell);
					});
				grid.Iterate_Faces(
					[&](const int axis, const VectorDi face) {
						const VectorDi cell_left = face;
				const VectorDi cell_right = face + VectorDi::Unit(1 - axis);
				initial_vel(axis, face) = (sol_host_big(cell_right) - sol_host_big(cell_left)) / grid.dx;
				if (axis == 1)
					initial_vel(axis, face) *= -1;

				if (face(axis) == 0 || face(axis) == grid.Counts()[axis])
					initial_vel(axis, face) = 0;
					});
				fluid.Init(cell_type, vol, face_fixed, initial_vel, true);
			}
		}
	};
}