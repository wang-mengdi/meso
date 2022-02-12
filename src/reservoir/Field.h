//////////////////////////////////////////////////////////////////////////
// Basic grid data (with data included)
// Copyright (c) (2022-), Bo Zhu, Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "Grid.h"

template<class T, int d>
class Field {
	Typedef_VectorD(d);
public:
	Grid<d, GridType::CELL> grid;
	Array<T> data;
	Field() {}
	Field(const Grid<d, GridType::CELL>& _grid, const real val = 0) :
		grid(_grid)
	{
		data.resize(grid.Size());
	}
};