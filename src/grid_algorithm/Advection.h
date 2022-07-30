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
	__global__ void RK1_Cell_Kernel(const real dt, const Grid<d> grid, T* result_data, const T* origin_data,
		const Grid<d> gv0, const T* v0, const Grid<d> gv1, const T* v1, const Grid<d> gv2, const T* v2) {
		Typedef_VectorD(d);
		Vector<int, d> cell = GPUFunc::Thread_Coord<d>(blockIdx, threadIdx);
		VectorD pos0 = grid.Position(cell);
		Vector<T, d> vel0 = Intp::Face_Vector(gv0, v0, gv1, v1, gv2, v2, pos0);
		VectorD back_pos = pos0 - vel0 * dt;
		result_data[grid.Index(cell)] = Intp::Value(grid, origin_data, back_pos);
	}

	template<class T, int d, class Intp>
	__global__ void RK2_Cell_Kernel(const real dt, const Grid<d> grid, T* result_data, const T* origin_data,
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

	template<class T, int d, class Intp>
	T RK2_Cell_Host(const real dt, const Field<T, d>& origin_field, const FaceField<T, d>& velocity, const Vector<int, d> cell) {
		Typedef_VectorD(d);
		VectorD pos0 = origin_field.grid.Position(cell);
		Vector<T, d> vel0 = Intp::Face_Vector(velocity, pos0);
		VectorD pos1 = pos0 - vel0 * 0.5 * dt;
		Vector<T, d> vel1 = Intp::Face_Vector(velocity, pos1);
		VectorD back_pos = pos0 - vel1 * dt;
		return Intp::Value(origin_field, back_pos);
	}

	template<class T, int d, class Intp>
	__global__ void RK3_Cell_Kernel(const real dt, const Grid<d> grid, T* result_data, const T* origin_data,
		const Grid<d> gv0, const T* v0, const Grid<d> gv1, const T* v1, const Grid<d> gv2, const T* v2) {
		Typedef_VectorD(d);
		Vector<int, d> cell = GPUFunc::Thread_Coord<d>(blockIdx, threadIdx);
		real c1 = 1.0 / 6.0 * dt, c2 = 4.0 / 6.0 * dt, c3 = 1.0 / 6.0 * dt;
		VectorD pos0 = grid.Position(cell);
		Vector<T, d> vel1 = Intp::Face_Vector(gv0, v0, gv1, v1, gv2, v2, pos0);
		VectorD pos1 = pos0 - vel1 * 0.5 * dt;
		Vector<T, d> vel2 = Intp::Face_Vector(gv0, v0, gv1, v1, gv2, v2, pos1);
		VectorD pos2 = pos0 + vel1 * dt - vel2 * 2.0 * dt;
		Vector<T, d> vel3 = Intp::Face_Vector(gv0, v0, gv1, v1, gv2, v2, pos2);
		VectorD back_pos = pos0 - c1 * vel1 - c2 * vel2 - c3 * vel3;
		result_data[grid.Index(cell)] = Intp::Value(grid, origin_data, back_pos);
	}

	template<class T, int d, class Intp>
	__global__ void RK4_Cell_Kernel(const real dt, const Grid<d> grid, T* result_data, const T* origin_data,
		const Grid<d> gv0, const T* v0, const Grid<d> gv1, const T* v1, const Grid<d> gv2, const T* v2) {
		Typedef_VectorD(d);
		Vector<int, d> cell = GPUFunc::Thread_Coord<d>(blockIdx, threadIdx);
		real c1 = 1.0 / 6.0 * dt, c2 = 1.0 / 3.0 * dt, c3 = 1.0 / 3.0 * dt, c4 = 1.0 / 6.0 * dt;
		VectorD pos0 = grid.Position(cell);
		Vector<T, d> vel1 = Intp::Face_Vector(gv0, v0, gv1, v1, gv2, v2, pos0);
		VectorD pos1 = pos0 - vel1 * 0.5 * dt;
		Vector<T, d> vel2 = Intp::Face_Vector(gv0, v0, gv1, v1, gv2, v2, pos1);
		VectorD pos2 = pos0 - vel2 * 0.5 * dt;
		Vector<T, d> vel3 = Intp::Face_Vector(gv0, v0, gv1, v1, gv2, v2, pos2);
		VectorD pos3 = pos0 - vel3 * dt;
		Vector<T, d> vel4 = Intp::Face_Vector(gv0, v0, gv1, v1, gv2, v2, pos3);
		VectorD back_pos = pos0 - c1 * vel1 - c2 * vel2 - c3 * vel3 - c4 * vel4;
		result_data[grid.Index(cell)] = Intp::Value(grid, origin_data, back_pos);
	}

	template<class Intp = IntpLinearPadding0, int order = 2>
	class Advection {
	public:
		template<class T, int d, DataHolder side>
		static void Advect(const real dt, Field<T, d, side>& advected_field, const Field<T, d, side>& origin_field, const FaceField<T, d, side>& velocity) {
			Assert(!std::is_same<Intp, IntpLinear>::value, "Advection can't use IntpLinear as interpolator");
			advected_field.Init(origin_field.grid);
			if constexpr (side == DEVICE) {
				const auto& vgrid = velocity.grid;
				Grid<d> vg0 = vgrid.Face_Grid(0), vg1 = vgrid.Face_Grid(1), vg2 = vgrid.Face_Grid(2);
				const T* v0 = velocity.Data_Ptr(0), * v1 = velocity.Data_Ptr(1), * v2 = velocity.Data_Ptr(2);
				if constexpr (order == 2)
					advected_field.grid.Exec_Kernel(
						&RK1_Cell_Kernel<T, d, Intp>, dt, advected_field.grid, advected_field.Data_Ptr(),
						origin_field.Data_Ptr(), vg0, v0, vg1, v1, vg2, v2);
				else if constexpr (order == 2)
					advected_field.grid.Exec_Kernel(
						&RK2_Cell_Kernel<T, d, Intp>, dt, advected_field.grid, advected_field.Data_Ptr(),
						origin_field.Data_Ptr(), vg0, v0, vg1, v1, vg2, v2);
				else if constexpr (order == 3)
					advected_field.grid.Exec_Kernel(
						&RK3_Cell_Kernel<T, d, Intp>, dt, advected_field.grid, advected_field.Data_Ptr(),
						origin_field.Data_Ptr(), vg0, v0, vg1, v1, vg2, v2);
				else if constexpr (order == 4)
					advected_field.grid.Exec_Kernel(
						&RK4_Cell_Kernel<T, d, Intp>, dt, advected_field.grid, advected_field.Data_Ptr(),
						origin_field.Data_Ptr(), vg0, v0, vg1, v1, vg2, v2);
				else
					Error("Invalid order for Advection.");
			}
			else {
				Error("Advection::Advect not implemented for HOST");
			}
		}

		//original_value and velocity can be the same
		//advected_result must be different
		//will allocate space for advected_result
		//don't handle boundary conditions
		template<class T, int d, DataHolder side>
		static void Advect(const real dt, FaceField<T, d, side>& advected_result, const FaceField<T, d, side>& original_value, const FaceField<T, d, side>& velocity) {
			advected_result.Init(original_value.grid);
			for (int axis = 0; axis < d; axis++) {
				Field<T, d, side> face_origin = original_value.Face_Reference(axis);
				Field<T, d, side> face_result = advected_result.Face_Reference(axis);
				Advect(dt, face_result, face_origin, velocity);
			}
		}
	};
}