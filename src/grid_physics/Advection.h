//////////////////////////////////////////////////////////////////////////
// Advection algorithms method
// Copyright (c) (2022-), Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once

#include "Field.h"
#include "FaceField.h"
#include "Interpolation.h"
#include "AuxFunc.h"

namespace Meso {

	template<class T, int d>
	__global__ static void Semi_Lagrangian_Cell(const real dt, const Grid<d> grid, T* result_data, const T* origin_data,
		const Grid<d> gv0, const T* v0, const Grid<d> gv1, const T* v1, const Grid<d> gv2, const T* v2) {
		Typedef_VectorD(d);
		Vector<int, d> cell = GPUFunc::Thread_Coord<d>(blockIdx, threadIdx);
		VectorD pos0 = grid.Position(cell);
		Vector<T, d> vel0 = IntpLinearPadding0::Face_Vector(gv0, v0, gv1, v1, gv2, v2, pos0);
		VectorD pos1 = pos0 - vel0 * 0.5 * dt;
		Vector<T, d> vel1 = IntpLinearPadding0::Face_Vector(gv0, v0, gv1, v1, gv2, v2, pos1);
		VectorD back_pos = pos0 - vel1 * dt;
		result_data[grid.Index(cell)] = IntpLinearPadding0::Value(grid, origin_data, back_pos);
	}

	template<int d>
	class SemiLagrangian {
		template<class T, DataHolder side>
		static void Advect(const real dt, Field<T, d, side>& advected_field, const FaceField<T, d, side>& velocity) {
			if constexpr (side == DEVICE) {
				Field<T, d, side> reference_field = advected_field;
				const auto& vgrid = velocity.grid;
				Grid<d> vg0 = vgrid.Face_Grid(0), vg1 = vgrid.Face_Grid(1), vg2 = vgrid.Face_Grid(2);
				const T* v0 = velocity.Data_Ptr(0), v1 = velocity.Data_Ptr(1), v2 = velocity.Data_Ptr(2);
				advected_field.Exec_Kernel(
					Semi_Lagrangian_Cell,
					dt,
					advected_field.grid,
					advected_field.Data_Ptr(),
					reference_field.Data_Ptr(),
					vg0, v0,
					vg1, v1,
					vg2, v2
				);
			}
			else {
				Error("Advection::Advect not implemented for HOST");
			}
		}

		template<class T, DataHolder side>
		static void Advect(FaceField<T, d, side>& advected_val, const FaceField<T, d, side>& velocity) {
			for (int axis = 0; axis < d; axis++) {
				Field<T, d, side> face_val(advected_val.grid.Face_Grid(axis), advected_val.face_data[axis]);
				Advect(face_val, velocity);
			}
		}
	};
}