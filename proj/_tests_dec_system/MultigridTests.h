//////////////////////////////////////////////////////////////////////////
// Multigrid operator and multigrid tests
// Copyright (c) (2022-), Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once

#include "Coarsener.h"
#include "Multigrid.h"
#include "PoissonTests.h"
#include "ConjugateGradient.h"
#include "Timer.h"

namespace Meso {

	void Test_Coarsener2(const Vector2i counts);
	void Test_Coarsener3(const Vector3i counts);

	template<class T, int d>
	void Test_Multigrid(const Vector<int, d> counts) {
		Check_Cuda_Memory("test begin");

		Grid<d> grid(counts);
		MaskedPoissonMapping<T, d> poisson = Random_Poisson_Mapping<T, d>(grid,1000);

		Field<T, d> b_host(grid);
		Random::Fill_Random_Array<T>(b_host.Data(), -5, 10);
		FieldDv<T, d> b_dev = b_host;
		FieldDv<T, d> x_dev(grid, 0);

		Check_Cuda_Memory("allocated poisson");

		FieldDv<T, d> res_dev(grid);
		poisson.Residual(res_dev.Data(), x_dev.Data(), b_dev.Data());
		Info("initial residual: {}", sqrt(ArrayFunc::Dot(res_dev.Data(), res_dev.Data())));

		VCycleMultigrid<T> solver;
		solver.Init_Poisson(poisson, 2, 2);

		solver.Apply(x_dev.Data(), b_dev.Data());
		poisson.Residual(res_dev.Data(), x_dev.Data(), b_dev.Data());
		Info("mg residual: {}", sqrt(ArrayFunc::Dot(res_dev.Data(), res_dev.Data())));

		//Info("run mg again:");
		//for (int k = 0; k < 10; k++) {
		//	solver.Apply(x_dev.Data(), b_dev.Data());
		//	poisson.Residual(res_dev.Data(), x_dev.Data(), b_dev.Data());
		//	Info("mg residual: {}", sqrt(ArrayFunc::Dot(res_dev.Data(), res_dev.Data())));
		//}
	}

	template<class T, int d>
	void Test_MGPCG(const Vector<int, d> counts) {
		Grid<d> grid(counts);
		MaskedPoissonMapping<T, d> poisson = Random_Poisson_Mapping<T, d>(grid);
		Field<T, d> b_host(grid);
		Random::Fill_Random_Array<T>(b_host.Data(), -5, 10);
		FieldDv<T, d> b_dev = b_host;
		FieldDv<T, d> x_dev(grid, 0);
		ConjugateGradient<T> MGPCG;
		VCycleMultigrid<T> precond;
		precond.Init_Poisson(poisson, 2, 2);
		//MGPCG.Init(&poisson, &precond, false);
		MGPCG.Init(&poisson, &precond, false, -1, 1e-6);
		//MGPCG.Init(&poisson, nullptr, true);
		int iters = 0;
		real res = 0;
		Timer timer;
		MGPCG.Solve(x_dev.Data(), b_dev.Data(), iters, res);
		//Info("MGPCG solved {} iters with relative_error={}", iters, res);
		if (iters < 100) {
			Pass("MGPCG test passed in {}s for counts={}, with {} iters and relative_error={}", timer.Lap_Time(), counts, iters, res);
		}
		else {
			Error("MGPCG test failed for counts={}, with {} iters and relative_error={}", counts, iters, res);
		}
	}

}