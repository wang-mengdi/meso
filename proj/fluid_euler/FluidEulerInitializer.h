//////////////////////////////////////////////////////////////////////////
// Initializer of a Fluid Euler System
// Copyright (c) (2022-), Bo Zhu, Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once

#include "FluidEuler.h"
#include "Json.h"

namespace Meso {
	template<int d>
	class FluidEulerInitializer {
	public:
		Typedef_VectorD(d);
		void Apply(json& j, FluidEuler<d>& fluid) {
			int test = Json::Value(j, "test", 0);
			switch (test) {
			case 0:Case_0(j, fluid); break;
			default:Assert(false, "test {} not exist", test); break;
			}
		}

		bool Cell_In_Boundary(const Grid<d> grid, const VectorDi cell, int axis, int side, int width) {
			if (side == 0) return cell[axis] < width;
			else if (side == 1) return cell[axis] >= grid.counts[axis] - width;
			return false;
		}

		bool Face_In_Boundary(const Grid<d> grid, int axis, const VectorDi face, int chk_axis, int side, int width) {
			int lbound = width;
			int rbound = (width == -1) ? grid.Face_Grid(axis).counts[chk_axis] + 1 : grid.counts[chk_axis] - width;
			if (side == 0) {
				if (axis == chk_axis) return face[chk_axis] <= lbound;
				else return face[chk_axis] < lbound;
			}
			else if (side == 1) {
				return face[chk_axis] >= rbound;
			}
			else return false;
		}

		//all these passed fields are not initialized before
		void Set_Boundary(const Grid<d> grid, const Eigen::Matrix<int, 3, 2> bc_width, const Eigen::Matrix<real, 3, 2> bc_val,
			Field<bool, d>& cell_fixed, FaceField<real, d>& vol, FaceField<bool, d> &face_fixed, FaceField<real, d>& boundary_vel) {
			cell_fixed.Init(grid);
			vol.Init(grid);
			face_fixed.Init(grid);
			boundary_vel.Init(grid);

			grid.Exec_Nodes(
				[&](const VectorDi cell) {
					cell_fixed(cell) = false;
					for (int axis = 0; axis < d; axis++) {
						for (int side = 0; side < 2; side++) {
							if (Cell_In_Boundary(grid, cell, axis, side, bc_width(axis,side))) {
								cell_fixed(cell) = true;
								return;
							}
						}
					}
				}
			);

			Info("grid counts: {}", grid.counts);

			grid.Exec_Faces(
				[&](const int axis, const VectorDi face) {
					vol(axis, face) = 1;
					face_fixed(axis, face) = false;
					boundary_vel(axis, face) = 0;
					int in_cnt = 0, _chk_axis, _side;
					for (int chk_axis = 0; chk_axis < d; chk_axis++) {
						for (int side = 0; side < 2; side++) {
							if (Face_In_Boundary(grid, axis, face, chk_axis, side, bc_width(chk_axis, side))) {
								in_cnt++;
								_chk_axis = chk_axis;
								_side = side;
							}
						}
					}
					if (in_cnt > 0) {
						if (in_cnt == 1 && axis == _chk_axis) {
							vol(axis, face) = 0;
							face_fixed(axis, face) = true;
							boundary_vel(axis, face) = bc_val(_chk_axis, _side);
							//Info("set face bc {},{} to {} with chk_axis {} and side {}", axis, face, boundary_vel(axis, face), _chk_axis, _side);
						}
						else {
							vol(axis, face) = 0;
							face_fixed(axis, face) = true;
							boundary_vel(axis, face) = 0;
						}
					}
				}
			);
		}

		//an empty field with inflow boundary condition v=1
		void Case_0(json& j, FluidEuler<d>& fluid) {
			int scale = Json::Value(j, "scale", 32);
			real dx = 1.0 / scale;
			VectorDi grid_size = scale * VectorFunc::Vi<d>(2, 1, 1);
			Grid<d> grid(grid_size, dx, VectorD::Zero(), MAC);
			Field<bool, d> fixed;
			FaceField<real, d> vol;
			FaceField<bool, d> face_fixed;
			FaceField<real, d> initial_vel;

			Eigen::Matrix<int, 3, 2> bc_width;
			Eigen::Matrix<real, 3, 2> bc_val;
			bc_width << 1, -1, -1, -1, -1, -1;
			bc_val << 1, 0, 0, 0, 0, 0;

			//Info("bc_width: \n{}\n, bc_val: \n{}\n", bc_width, bc_val);
			std::cout << "bc_width: " << bc_width << "\n";
			std::cout << "bc_width: " << bc_val << "\n";

			Set_Boundary(grid, bc_width, bc_val, fixed, vol, face_fixed, initial_vel);
			fluid.Init(fixed, vol, face_fixed, initial_vel);
			//ArrayFunc::Fill(fluid.velocity.Data(0), 1.0);
		}
	};
}