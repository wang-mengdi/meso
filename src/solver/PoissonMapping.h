//////////////////////////////////////////////////////////////////////////
// Poisson Linear Mapping
// Copyright (c) (2018-), Zangyueyang Xian, Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "LinearMapping.h"
#include "Field.h"
#include "LambdaHelper.h"

template<class T, int d>
class PoissonMapping : public LinearMapping<T> {
	Typedef_VectorD(d);
public:
	FieldDv<T, d> vol;
	FieldDv<bool, d> fixed;

	virtual int xDoF() const {}//number of cols

	virtual int yDoF() const {}//number of rows

	//input p, get Ap
	virtual void applyMapping(ArrayDv<T>& Ap, const ArrayDv<T>& p) {}

	void Init(const Grid<d, GridType::CELL>& grid, IFFunc<T, d> vol_func, CFunc<T, d> is_unknown_func) {
		int dof = grid.DoF();

	}
};