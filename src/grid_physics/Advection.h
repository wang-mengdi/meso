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
#include "BoundaryCondition.h"

namespace Meso {

	template<class T, int d, class Intp>
	__global__ static void Semi_Lagrangian_Cell(const real dt, const Grid<d> grid, T* result_data, const T* origin_data,
		const Grid<d> gv0, const T* v0, const Grid<d> gv1, const T* v1, const Grid<d> gv2, const T* v2) {
		Typedef_VectorD(d);
		Vector<int, d> cell = GPUFunc::Thread_Coord<d>(blockIdx, threadIdx);
		VectorD pos0 = grid.Position(cell);
		Vector<T, d> vel0 = Intp::Face_Vector(gv0, v0, gv1, v1, gv2, v2, pos0);
		VectorD pos1 = pos0 - vel0 * 0.5 * dt;
		Vector<T, d> vel1 = Intp::Face_Vector(gv0, v0, gv1, v1, gv2, v2, pos1);
		VectorD back_pos = pos0 - vel1 * dt;
		result_data[grid.Index(cell)] = Intp::Value(grid, origin_data, back_pos);
	}

	template<class T, int d>
	__global__ static void Inver_Flow_Map_Cell(const real dt, const Grid<d> grid, Vector<T,d>* inverse_flow_map,
		const Grid<d> gv0, const T* v0, const Grid<d> gv1, const T* v1, const Grid<d> gv2, const T* v2) {
		Typedef_VectorD(d);
		Vector<int, d> cell = GPUFunc::Thread_Coord<d>(blockIdx, threadIdx);
		VectorD pos0 = grid.Position(cell);
		Vector<T, d> vel0 = IntpLinearPadding0::Face_Vector(gv0, v0, gv1, v1, gv2, v2, pos0);
		VectorD pos1 = pos0 - vel0 * 0.5 * dt;
		Vector<T, d> vel1 = IntpLinearPadding0::Face_Vector(gv0, v0, gv1, v1, gv2, v2, pos1);
		VectorD back_pos = pos0 - vel1 * dt;
		inverse_flow_map[grid.Index(cell)] = back_pos;
	}

	template<class Intp = IntpLinearPadding0>
	class SemiLagrangian {
	public:
		template<class T, int d, DataHolder side>
		static void Advect(const real dt, Field<T, d, side>& advected_field, const Field<T, d, side>& origin_field, const FaceField<T, d, side>& velocity) {
			advected_field.Init(origin_field.grid);
			if constexpr (side == DEVICE) {
				const auto& vgrid = velocity.grid;
				Grid<d> vg0 = vgrid.Face_Grid(0), vg1 = vgrid.Face_Grid(1), vg2 = vgrid.Face_Grid(2);
				const T* v0 = velocity.Data_Ptr(0), * v1 = velocity.Data_Ptr(1), * v2 = velocity.Data_Ptr(2);
				advected_field.grid.Exec_Kernel(
					&Semi_Lagrangian_Cell<T, d, Intp>,
					dt,
					advected_field.grid,
					advected_field.Data_Ptr(),
					origin_field.Data_Ptr(),
					vg0, v0,
					vg1, v1,
					vg2, v2
				);
			}
			else {
				Error("Advection::Advect not implemented for HOST");
			}
		}

		template<class T, int d, DataHolder side>
		static void Inverse_Flow_Map(const real dt, Field<Vector<real, d>, d, side>& inverse_flow_map, const FaceField<T, d, side>& velocity) {
			if constexpr (side == DEVICE) {
				const auto& vgrid = velocity.grid;
				Grid<d> vg0 = vgrid.Face_Grid(0), vg1 = vgrid.Face_Grid(1), vg2 = vgrid.Face_Grid(2);
				const T* v0 = velocity.Data_Ptr(0), * v1 = velocity.Data_Ptr(1), * v2 = velocity.Data_Ptr(2);
				inverse_flow_map.grid.Exec_Kernel(
					&Inver_Flow_Map_Cell<T, d>,
					dt,
					inverse_flow_map.grid,
					inverse_flow_map.Data_Ptr(),
					vg0, v0,
					vg1, v1,
					vg2, v2
				);
			}
			else {
				Error("Advection::Inverse_Flow_Map not implemented for HOST");
			}
		}

		////advected_val may be the same as velocity
		//template<class T, DataHolder side>
		//static void Advect(const real dt, FaceField<T, d, side>& advected_val, const FaceField<T, d, side>& velocity, const BoundaryCondition<FaceField<T, d, side>>& bc) {
		//	FaceField<T, d, side> advection_result(advected_val.grid);
		//	for (int axis = 0; axis < d; axis++) {
		//		const auto face_grid = advected_val.grid.Face_Grid(axis);
		//		Field<T, d, side> face_origin(face_grid, advected_val.face_data[axis]);
		//		Field<T, d, side> face_result(face_grid, advection_result.face_data[axis]);
		//		Advect(dt, face_result, face_origin, velocity);
		//	}
		//	bc.Copy_UnMasked(advected_val, advection_result);
		//}

		//original_value and velocity can be the same
		//advected_result must be different
		//will allocate space for advected_result
		//don't handle boundary conditions
		template<class T, int d, DataHolder side>
		static void Advect(const real dt, FaceField<T, d, side>& advected_result, const FaceField<T, d, side>& original_value, const FaceField<T, d, side>& velocity) {
			advected_result.Init(original_value.grid);
			for (int axis = 0; axis < d; axis++) {
				const auto face_grid = original_value.grid.Face_Grid(axis);
				Field<T, d, side> face_origin(face_grid, original_value.face_data[axis]);
				Field<T, d, side> face_result(face_grid, advected_result.face_data[axis]);
				Advect(dt, face_result, face_origin, velocity);
			}
		}
	};
}