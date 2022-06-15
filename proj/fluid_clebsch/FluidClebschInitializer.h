//////////////////////////////////////////////////////////////////////////
// Initializer of a Clebsch Fluid System
// Copyright (c) (2022-), Bo Zhu, Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once

#include "FluidClebsch.h"
#include "ImplicitGeometry.h"
#include "Json.h"
#include "GridEulerFunc.h"
#include "AuxFunc.h"

namespace Meso {
	template<int d>
	class  FluidClebschInitializer{
	public:
		Typedef_VectorD(d);

		//define the boundary conditions of the eulerian system
		Field<bool, d> fixed;
		Field<Vector2C, d> initial_wave_func;
		FaceField<real, d> vol;

		void Apply(json& j, FluidClebsch<d>& fluid) {
			int test = Json::Value(j, "test", 0);
			switch (test) {
			case 0:Case_0(j, fluid); break;
			default:Assert(false, "test {} not exist", test); break;
			}
		}

		void Set_Boundary(const Grid<d> grid, const Eigen::Matrix<int, 3, 2> bc_width, const Eigen::Matrix<real, 3, 2> bc_val,
			Field<bool, d>& cell_fixed, FaceField<real, d>& vol) {
			cell_fixed.Init(grid);
			vol.Init(grid);

			grid.Exec_Nodes(
				[&](const Vector<int, d> cell) {
					cell_fixed(cell) = false;
					for (int axis = 0; axis < d; axis++) {
						for (int side = 0; side < 2; side++) {
							if (GridEulerFunc::Cell_In_Boundary<d>(grid, cell, axis, side, bc_width(axis, side))) {
								cell_fixed(cell) = true;
								return;
							}
						}
					}
				}
			);

			grid.Exec_Faces(
				[&](const int axis, const Vector<int, d> face) {
					vol(axis, face) = 1;
					int in_cnt = 0, _chk_axis, _side;
					for (int chk_axis = 0; chk_axis < d; chk_axis++) {
						for (int side = 0; side < 2; side++) {
							if (GridEulerFunc::Face_In_Boundary<d>(grid, axis, face, chk_axis, side, bc_width(chk_axis, side))) {
								in_cnt++;
								_chk_axis = chk_axis;
								_side = side;
							}
						}
					}
					if (in_cnt > 0) {
						if (in_cnt == 1 && axis == _chk_axis) {
							vol(axis, face) = 0;
						}
						else {
							vol(axis, face) = 0;
						}
					}
				}
			);
		}

		void Initialize_Vortex_Ring(const Grid<d> grid, Field<Vector2C, d>& initial_wave_func, const VectorD velocity) {
			initial_wave_func.Init(grid);
			const VectorD c = grid.Center();
			const real r = (real)0.1;
			grid.Exec_Nodes(
				[&](const Vector<int, d> cell) {
					const VectorD pos = grid.Position(cell);
					Vector2C psi(thrust::complex<real>{1.0 / std::sqrt(1.01), 0.}, thrust::complex<real>{0.1 / std::sqrt(1.01), 0.});
					const real phase = velocity.dot(pos);
					psi[0] *= thrust::complex<real>{std::exp(1i * phase)};
					psi[1] *= thrust::complex<real>{std::exp(1i * phase)};
					const real rx = (pos[0] - c[0]) / r;
					const real r2 = (pos - c).squaredNorm() / std::pow(r, 2);
					const real fR = std::exp(-std::pow(r2 / (real)9, (real)4));
					const thrust::complex<real> q(2. * rx * fR / (r2 + 1), (r2 + 1. - 2. * fR) / (r2 + 1));
					psi[0] *= q;
					initial_wave_func(cell) = psi;
				}
			);
		}

		//initialized vortex ring projection test
		void Case_0(json& j, FluidClebsch<d>& fluid) {
			int scale = Json::Value(j, "scale", 32);
			real dx = 1.0 / scale;
			VectorDi grid_size = scale * MathFunc::Vi<d>(1, 1, 1);
			Grid<d> grid(grid_size, dx, VectorD::Zero(), MAC);

			Eigen::Matrix<int, 3, 2> bc_width;
			Eigen::Matrix<real, 3, 2> bc_val;
			bc_width << 0, -1, 0, 0, 0, 0;
			bc_val << 0, 0, 0, 0, 0, 0;

			Set_Boundary(grid, bc_width, bc_val, fixed, vol);
			VectorD background_velocity = 0. * VectorD::Unit(1);
			Initialize_Vortex_Ring(grid, initial_wave_func, background_velocity);
			real h_bar = Json::Value(j, "h_bar", 0.01);
			fluid.Init(h_bar, fixed, initial_wave_func, vol);
		}
	};
}