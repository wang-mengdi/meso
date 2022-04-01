//////////////////////////////////////////////////////////////////////////
// Interpolation on grid
// Copyright (c) (2022-), Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once

#include "Grid.h"

namespace Meso {
	namespace Interpolation {
		template<class T, int d, GridType gtype>
		T Linear_Intp(const Grid<d, gtype> grid, const T* data, const Vector<int, d> coord, const Vector<real, d> frac) {
			Typedef_VectorD(d);
			//if constexpr (d == 1) {
			//	auto a0l = 1.0 - frac[0], a0r = frac[0];
			//	int i = cell[0];
			//	int indexl = grid.Index(i);
			//	int indexr = grid.Index(i + 1);
			//	return
			//		a0l * data[indexl] +
			//		a0r * data[indexr];
			//}
			if constexpr (d == 2) {
				auto a0l = 1.0 - frac[0], a0r = frac[0];
				auto a1l = 1.0 - frac[1], a1r = frac[1];
				int i = coord[0], j = coord[1];
				int indexll = grid.Index(i, j);
				int indexlr = grid.Index(i, j + 1);
				int indexrl = grid.Index(i + 1, j);
				int indexrr = grid.Index(i + 1, j + 1);
				return
					a0l * a1l * data[indexll] +
					a0l * a1r * data[indexlr] +
					a0r * a1l * data[indexrl] +
					a0r * a1r * data[indexrr];
			}
			else if constexpr (d == 3) {
				auto a0l = 1.0 - frac[0], a0r = frac[0];
				auto a1l = 1.0 - frac[1], a1r = frac[1];
				auto a2l = 1.0 - frac[2], a2r = frac[2];
				int i = coord[0], j = coord[1], k = coord[2];
				int indexlll = grid.Index(i, j, k);
				int indexllr = grid.Index(i, j, k + 1);
				int indexlrl = grid.Index(i, j + 1, k);
				int indexlrr = grid.Index(i, j + 1, k + 1);
				int indexrll = grid.Index(i + 1, j, k);
				int indexrlr = grid.Index(i + 1, j, k + 1);
				int indexrrl = grid.Index(i + 1, j + 1, k);
				int indexrrr = grid.Index(i + 1, j + 1, k + 1);
				return
					a0l * a1l * a2l * data[indexlll] +
					a0l * a1l * a2r * data[indexllr] +
					a0l * a1r * a2l * data[indexlrl] +
					a0l * a1r * a2r * data[indexlrr] +
					a0r * a1l * a2l * data[indexrll] +
					a0r * a1l * a2r * data[indexrlr] +
					a0r * a1r * a2l * data[indexrrl] +
					a0r * a1r * a2r * data[indexrrr];
			}
			else Assert("Interpolation:Linear_Intp error: dimension must be 2 or 3");
		}
	}
}