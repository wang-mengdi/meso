//////////////////////////////////////////////////////////////////////////
// Common auxillary functions
// Copyright (c) (2022-), Bo Zhu, Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "Common.h"

namespace VectorFunc {
	////create vectors with compatible dimensions
	template<int d> Vector<real, d> V(const real x = (real)0, const real y = (real)0, const real z = (real)0);
	template<int d> Vector<int, d> Vi(const int x = 0, const int y = 0, const int z = 0, const int w = 0);
	template<int d> Vector<real, d> V(const Vector2 v2);
	template<int d> Vector<real, d> V(const Vector3 v2);

	////Round up vector to a multiple of bn
	template<int d> Vector<int, d> Round_Up_To_Align(Vector<int, d> v, int bn) {
		for (int i = 0; i < d; i++) {
			//round to the nearest multiple of block_size
			v[i] = ((v[i] + bn - 1) / bn) * bn;
		}
		return v;
	}
}