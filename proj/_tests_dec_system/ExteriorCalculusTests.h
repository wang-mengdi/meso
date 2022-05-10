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
			T v0 = (grid.Valid(cell0) ? field_host(cell0) : 0);
			T v1 = (grid.Valid(cell1) ? field_host(cell1) : 0);
			return v1 - v0;
		}
	);
	for (int axis = 0; axis < d; axis++) {
		if (!ArrayFunc::Is_Approx<T>(facefield_naive.Data(axis), facefield_ext_host.Data(axis))) {
			Error("Test_Exterior_Derivative_Cell failed for counts={}", counts);
			return;
		}
	}
	Pass("Test_Exterior_Derivative_Cell passed for counts={}", counts);
}

template<class T, int d>
void Test_Exterior_Derivative_Face(const Vector<int, d> counts) {
	Typedef_VectorD(d);
	Grid<d> grid(counts);
	FaceField<T, d> F_host(grid);
	for (int axis = 0; axis < d; axis++) Random::Fill_Random_Array<T>(F_host.Data(axis), -3, 10);
	//F_host.Calc_Faces(
	//	[&](const int axis, const VectorDi face) {
	//		if (axis == 0) return 1;
	//		else return 0;
	//	}
	//);
	FaceFieldDv<T, d> F_dev = F_host;
	FieldDv<T, d> C_ext_dev;
	Exterior_Derivative(C_ext_dev, F_dev);
	Field<T, d> C_ext_host = C_ext_dev;
	Field<T, d> C_naive(grid);
	C_naive.Calc_Cells(
		[&](const VectorDi cell) {
			T div = 0;
			for (int axis = 0; axis < d; axis++) {
				VectorDi face0 = cell, face1 = cell + VectorDi::Unit(axis);
				div += F_host(axis, face1) - F_host(axis, face0);
			}
			return div;
		}
	);
	if (ArrayFunc::Is_Approx<T>(C_naive.Data(), C_ext_host.Data())) {
		Pass("Test_Exterior_Derivative_Face passed for counts={}", counts);
	}
	else {
		Error("Test_Exterior_Derivative_Face passed for counts={}", counts);
	}
}