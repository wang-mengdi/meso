//////////////////////////////////////////////////////////////////////////
// Advection algorithms method
// Copyright (c) (2022-), Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once

#include "Field.h"
#include "FaceField.h"

namespace Meso {
	template<int d>
	class SemiLagrangian {
		template<class T, DataHolder side, GridType gtype>
		__global__ static void Semi_Lagrangian_Cell(Field<T, d, side, gtype>& advected_field, const Field<T, d, side, gtype>& advected_field, const FaceField<T, d, side>& velocity) {
			
		}

		template<class T, DataHolder side, GridType gtype>
		static void Advect(Field<T, d, side, gtype>& advected_field, const FaceField<T, d, side>& velocity) {
			if constexpr (side == DEVICE) {
				Field<T, d, side, gtype> reference_field = advected_field;
				advected_field.Exec_Kernel(
				);
			}
			else {
				Error("Advection::Advect not implemented for HOST");
			}
		}

		template<class T, DataHolder side>
		static void Advect(FaceField<T, d, side>& advected_val, const FaceField<T, d, side>& velocity) {
			for (int axis = 0; axis < d; axis++) {
				Field<T, d, side, CORNER> face_val(advected_val.grid.Face_Grid(axis), advected_val.face_data[axis]);
				Advect(face_val, velocity);
			}
		}
	};
}