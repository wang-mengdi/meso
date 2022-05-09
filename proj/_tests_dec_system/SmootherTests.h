//////////////////////////////////////////////////////////////////////////
// Smoother Tests
// Copyright (c) (2022-), Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "GridGSSmoother.h"
#include "DampedJacobiSmoother.h"
#include "PoissonTests.h"

namespace Meso {

	template<class T, int d>
	void Test_GridGSSmoother(const Vector<int, d> counts) {
		Grid<d> grid(counts);
		MaskedPoissonMapping<T, d> poisson = Random_Poisson_Mapping<T, d>(grid);
		Field<T, d> b_host(grid);
		Random::Fill_Random_Array<T>(b_host.Data(), -5, 10);
		FieldDv<T, d> b_dev = b_host;
		FieldDv<T, d> x_dev(grid, 0);
		FieldDv<T, d> res_dev(grid);
		poisson.Residual(res_dev.Data(), x_dev.Data(), b_dev.Data());
		Info("initial residual: {}", sqrt(ArrayFunc::Dot(res_dev.Data(), res_dev.Data())));

		for (int iter_num = 0; iter_num < 10; iter_num++) {
			GridGSSmoother<T, d> smoother(poisson, iter_num);
			smoother.Apply(x_dev.Data(), b_dev.Data());
			poisson.Residual(res_dev.Data(), x_dev.Data(), b_dev.Data());
			Info("iter_num={} residual={}", iter_num, sqrt(ArrayFunc::Dot(res_dev.Data(), res_dev.Data())));
		}

	}

	template<class T, int d>
	void Test_DampedJacobiSmoother(const Vector<int, d> counts) {
		Grid<d> grid(counts);
		MaskedPoissonMapping<T, d> poisson = Random_Poisson_Mapping<T, d>(grid);
		Field<T, d> b_host(grid);
		Random::Fill_Random_Array<T>(b_host.Data(), -5, 10);
		FieldDv<T, d> b_dev = b_host;
		FieldDv<T, d> x_dev(grid, 0);
		FieldDv<T, d> res_dev(grid);
		poisson.Residual(res_dev.Data(), x_dev.Data(), b_dev.Data());
		Info("initial residual: {}", sqrt(ArrayFunc::Dot(res_dev.Data(), res_dev.Data())));
		ArrayDv<T> poisson_diag; Poisson_Diagonal(poisson_diag, poisson);

		for (int iter_num = 0; iter_num < 10; iter_num++) {
			DampedJacobiSmoother<T> smoother(poisson, poisson_diag, iter_num, 2.0 / 3.0);
			smoother.Apply(x_dev.Data(), b_dev.Data());
			poisson.Residual(res_dev.Data(), x_dev.Data(), b_dev.Data());
			Info("iter_num={} residual={}", iter_num, sqrt(ArrayFunc::Dot(res_dev.Data(), res_dev.Data())));
		}

	}

}