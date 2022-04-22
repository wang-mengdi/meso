//////////////////////////////////////////////////////////////////////////
// Test Differential Exterior Calculus
// Copyright (c) (2022-), Zangyueyang Xian, Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once

#include "DifferentialExteriorCalculus.h"
#include "Random.h"
using namespace Meso;

//test: facefield=d(field)
template<class T, int d>
void Test_Exterior_Derivative_Cell(const Vector<int, d> counts) {
	Typedef_VectorD(d);
	Grid<d> grid(counts);
	Field<T, d> field_host(grid);
	Random::Fill_Random_Array<T>(field_host.Data(), -3, 10);
	FieldDv<T, d> field_dev = field_host;
	FaceFieldDv<T, d> facefield_ext_dev;
	Exterior_Derivative(facefield_ext_dev, field_dev);
	FaceField<T, d> facefield_ext_host = facefield_ext_dev;
	FaceField<T, d> facefield_naive(grid);
	facefield_naive.Calc_Faces(
		[&](const int axis, const VectorDi face) {
			VectorDi cell0 = face - VectorDi::Unit(axis);
			VectorDi cell1 = face;
			real v0 = (grid.Valid(cell0) ? field_host(cell0) : 0);
			real v1 = (grid.Valid(cell1) ? field_host(cell1) : 0);
			return v1 - v0;
		}
	);
	for (int axis = 0; axis < d; axis++) {
		if (!ArrayFunc::IsApprox<T>(facefield_naive.Data(axis), facefield_ext_host.Data(axis))) {
			Error("Test_Exterior_Derivative_Cell failed for counts={}", counts);
			return;
		}
	}
	Pass("Test_Exterior_Derivative_Cell passed for counts={}", counts);
}