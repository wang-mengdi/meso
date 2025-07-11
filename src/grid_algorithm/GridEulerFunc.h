//////////////////////////////////////////////////////////////////////////
// Grid Functions
// Copyright (c) (2022-), Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once

#include "Grid.h"
#include "Field.h"
#include "FaceField.h"

namespace Meso {
	namespace GridEulerFunc {
		template<class T, int d, DataHolder side>
		T Linf_Norm(Field<T, d, side>& F) {
			F.Set_Padding_To(0);
			return ArrayFunc::Max_Abs<T, side>(F.Data());
		}

		template<class T, int d, DataHolder side>
		T Linf_Norm(Field<Vector<T,d>, d, side>& F) {
			F.Set_Padding_To(Vector<T,d>::Zero());
			return ArrayFunc::Largest_Norm<Vector<T, d>>(F.Data());
		}

		template<class T, int d, DataHolder side>
		T Linf_Norm(FaceField<T, d, side>& F) {
			T axis_max = (T)0;
			for (int axis = 0; axis < d; axis++) {
				Field<T, d, side>  field = F.Face_Reference(axis);
				T axis_norm = Linf_Norm(field);
				axis_max = axis_norm > axis_max ? axis_norm : axis_max;
			}
			return axis_max;
		}

		template<int d>
		bool Cell_In_Boundary(const Grid<d> grid, const Vector<int, d> cell, int axis, int side, int width) {
			if (side == 0) return cell[axis] < width;
			else if (side == 1) return cell[axis] >= grid.Counts()[axis] - width;
			return false;
		}

		template<int d>
		bool Face_In_Boundary(const Grid<d> grid, int axis, const Vector<int, d> face, int chk_axis, int side, int width) {
			int lbound = width;
			int rbound = (width == -1) ? grid.Face_Grid(axis).Counts()[chk_axis] + 1 : grid.Counts()[chk_axis] - width;
			if (side == 0) {
				if (axis == chk_axis) return face[chk_axis] <= lbound;
				else return face[chk_axis] < lbound;
			}
			else if (side == 1) {
				return face[chk_axis] >= rbound;
			}
			else return false;
		}

		template<class T, int d>
		void Set_Boundary_Cells(Field<T, d>& F, const Eigen::Matrix<int, 3, 2> bc_width, const T val) {
			F.Exec_Nodes(
				[&](const Vector<int, d> cell) {
					for (int axis = 0; axis < d; axis++) {
						for (int side = 0; side < 2; side++) {
							if (Cell_In_Boundary<d>(F.grid, cell, axis, side, bc_width(axis, side))) {
								F(cell) = val;
								return;
							}
						}
					}
				}
			);
		}

		//all these passed fields are not initialized before
		template<int d>
		void Set_Boundary(const Grid<d> grid, const Eigen::Matrix<int, 3, 2> bc_width, const Eigen::Matrix<real, 3, 2> bc_val,
			Field<unsigned char, d>& cell_type, FaceField<real, d>& vol, FaceField<bool, d>& face_fixed, FaceField<real, d>& boundary_vel) {
			
			cell_type.Init(grid, 0);
			Set_Boundary_Cells<unsigned char>(cell_type, bc_width, 2);

			vol.Init(grid);
			face_fixed.Init(grid);
			boundary_vel.Init(grid);
			grid.Exec_Faces(
				[&](const int axis, const Vector<int, d> face) {
					vol(axis, face) = 1;
					face_fixed(axis, face) = false;
					boundary_vel(axis, face) = 0;
					int in_cnt = 0, _chk_axis, _side;
					for (int chk_axis = 0; chk_axis < d; chk_axis++) {
						for (int side = 0; side < 2; side++) {
							if (Face_In_Boundary<d>(grid, axis, face, chk_axis, side, bc_width(chk_axis, side))) {
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
	}
}
