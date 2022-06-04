//////////////////////////////////////////////////////////////////////////
// Sparse solver tests
// Copyright (c) (2022-), Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "SparseDirectSolver.h"
#include "PoissonFunc.h"
#include "PoissonTests.h"

namespace Meso {
	template<class T, int d>
	void Test_Cholesky_Sparse_Solve(const Vector<int, d> counts) {
		Grid<d> grid(counts);
		MaskedPoissonMapping<T, d> poisson_mapping = Random_Poisson_Mapping<T, d>(grid, 1000);

		//Array<T> p_host = Random::Random_Array<T>(sparse_mapping.YDoF(), 0, 10);
		//ArrayDv<T> p_dev = p_host;
		//ArrayDv<T> Ap_poisson(poisson_mapping.XDoF()), Ap_sparse(sparse_mapping.XDoF());
		//poisson_mapping.Apply(Ap_poisson, p_dev);
		//sparse_mapping.Apply(Ap_sparse, p_dev);
		//Assert(ArrayFunc::Is_Approx<T>(Ap_poisson, Ap_sparse), "sparse mapping not correct");
		//Info("p_dev: {}", p_dev);
		//Info("Ap_sparse: {}", Ap_sparse);
		//Info("Ap_poisson: {}", Ap_poisson);

		CholeskySparseSolver<T> solver(SparseMatrix_From_PoissonLike(grid, poisson_mapping));
		Array<T> b_host = Random::Random_Array<T>(poisson_mapping.YDoF(), 0, 10);
		ArrayDv<T> b_dev = b_host;
		ArrayDv<T> x_dev(poisson_mapping.YDoF()), res;
		solver.Apply(x_dev, b_dev);
		poisson_mapping.Residual(res, x_dev, b_dev);
		T b2 = ArrayFunc::Dot(b_dev, b_dev);
		T res2 = ArrayFunc::Dot(res, res);
		if (res2 / b2 < 1e-9) {
			Pass("Test_Cholesky_Sparse_Solve passed for {}", counts);
		}
		else {
			Error("Test_Cholesky_Sparse_Solve failed for {}", counts);
		}
	}
}