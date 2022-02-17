//////////////////////////////////////////////////////////////////////////
// Face data on MacGrid
// Copyright (c) (2022-), Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "Grid.h"

namespace Meso {

	template<class T, int d, DataHolder side = DataHolder::HOST>
	class FaceField {
		Typedef_VectorD(d);
	public:
		Grid<d, GridType::CELL> grid;
		std::array<Array<T, side>, d> face_data;
		FaceField() {}
		FaceField(const Grid<d, GridType::CELL>& _grid)
		{
			Init(_grid);
		}
		void Init(const Grid<d, GridType::CELL>& _grid) {
			grid = _grid;
			for (int axis = 0; axis < d; axis++) face_data[axis].resize(grid.DoF());
		}
	};

	template<class T, int d> using FaceFieldDv = FaceField<T, d, DEVICE>;

}