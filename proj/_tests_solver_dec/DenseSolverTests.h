//////////////////////////////////////////////////////////////////////////
// Dense matrix mapping tests
// Copyright (c) (2022-), Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "LUDirectSolver.h"
#include "PoissonTests.h"
using namespace Meso;

template<class T, int d>
void Test_LU_Dense_Solver(const Vector<int, d> counts) {
	Grid<d> grid(counts);
	MaskedPoissonMapping<T, d> poisson_mapping = Random_Poisson_Mapping<T, d>(grid);
	DenseMatrixMapping<T> dense_mapping;
	DenseMatrixMapping_From_PoissonLike(dense_mapping, grid, poisson_mapping);
	//dense_mapping.Init_PoissonLike(grid, poisson_mapping);
	LUDenseSolver<T> solver;
	solver.Init(dense_mapping);
	int n = dense_mapping.YDoF();
	Array<T> b_host = Random::Random_Array<T>(n);
	ArrayDv<T> b_dev(n); ArrayFunc::Copy(b_dev, b_host);
	ArrayDv<T> x_dev(n);
	ArrayDv<T> res(n);
	solver.Apply(x_dev, b_dev);
	poisson_mapping.Residual(res, x_dev, b_dev);
	T b2 = ArrayFunc::Dot(b_dev, b_dev);
	T res2 = ArrayFunc::Dot(res, res);
	if (res2 / b2 < 1e-9) {
		Pass("Test_LU_Dense_Solver passed for {}", counts);
	}
	else {
		Error("Test_LU_Dense_Solver failed for {}", counts);
	}

	//Info("run dense again:");
	//for (int i = 0; i < 100000000; i++) {
	//	solver.Apply(x_dev, b_dev);
	//	poisson_mapping.Residual(res, x_dev, b_dev);
	//	Info("iter {} res norm2: {}", i, ArrayFunc::Dot(res, res));
	//}
}