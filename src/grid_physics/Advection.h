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
	class Advection {
		template<class T, DataHolder side, GridType gtype>
		static void Advect(Field<T, d, side, gtype>& advected_field, const FaceField<T, d, side>& velocity) {

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