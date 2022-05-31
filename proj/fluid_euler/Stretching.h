#include "Field.h"
#include "FaceField.h"
#include "Interpolation.h"
#include "AuxFunc.h"
#include "BoundaryCondition.h"
#pragma once

namespace Meso {
	template<class T, int d>
	__global__ static void Covector_Stretching_Cell(const Grid<d> grid, T* result_data, const Vector<T, d>* inverse_flow_map_grad, 
		const Grid<d> gv0, const T* v0, const Grid<d> gv1, const T* v1, const Grid<d> gv2, const T* v2) {
		Typedef_VectorD(d);
		Vector<int, d> cell = GPUFunc::Thread_Coord<d>(blockIdx, threadIdx);
		VectorD pos = grid.Position(cell);
		Vector<T, d> grad = inverse_flow_map_grad[grid.Index(cell)];
		Vector<T, d> vel = IntpLinearPadding0::Face_Vector(gv0, v0, gv1, v1, gv2, v2, pos);
		result_data[grid.Index(cell)] = grad.dot(vel);
	}

	namespace Stretching {
		template<class T, int d, DataHolder side>
		void Covector_Stretching(Field<T, d, side>& stretched_field, const Field<Vector<T, d>, d, side>& inverse_flow_map_grad_field, const FaceField<T, d, side>& velocity) {
			if constexpr (side == DEVICE) {
				const auto& vgrid = velocity.grid;
				Grid<d> vg0 = vgrid.Face_Grid(0), vg1 = vgrid.Face_Grid(1), vg2 = vgrid.Face_Grid(2);
				const T* v0 = velocity.Data_Ptr(0), * v1 = velocity.Data_Ptr(1), * v2 = velocity.Data_Ptr(2);
				stretched_field.grid.Exec_Kernel(
					&Covector_Stretching_Cell<T, d>,
					inverse_flow_map_grad_field.grid,
					stretched_field.Data_Ptr(),
					inverse_flow_map_grad_field.Data_Ptr(),
					vg0, v0,
					vg1, v1,
					vg2, v2
				);
			}
			else {
				Error("Stretching not implemented for HOST");
			}
			
		}

		template<class T, int d, DataHolder side>
		void Covector_Stretching(FaceField<T, d, side>& stretched_result, const FaceField<Vector<T,d>, d, side>& inverse_flow_map_grad, const FaceField<T, d, side>& velocity) {
			for (int axis = 0; axis < d; axis++) {
				const auto face_grid = inverse_flow_map_grad.grid.Face_Grid(axis);
				Field<Vector<T, d>, d, side> inverse_flow_map_grad_field(face_grid, inverse_flow_map_grad.face_data[axis]);
				Field<T, d, side> face_result(face_grid, stretched_result.face_data[axis]);
				Covector_Stretching(face_result, inverse_flow_map_grad_field, velocity);
			}
		}
	}
}