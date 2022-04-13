//////////////////////////////////////////////////////////////////////////
// Multigrid tests
// Copyright (c) (2022-), Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once

#include "Multigrid.h"
#include "PoissonTests.h"
#include "ConjugateGradient.h"
using namespace Meso;

template<class T, int d>
void Test_Multigrid(const Vector<int, d> counts) {
	Check_Cuda_Memory("test begin");

	Grid<d> grid(counts);
	PoissonMapping<T, d> poisson = Random_Poisson_Mapping<T, d>(grid);

	Field<T, d> b_host(grid);
	Random::Fill_Random_Array<T>(b_host.data, -5, 10);
	FieldDv<T, d> b_dev = b_host;
	FieldDv<T, d> x_dev(grid, 0);

	Check_Cuda_Memory("allocated poisson");

	FieldDv<T, d> res_dev(grid);
	poisson.Residual(res_dev.data, x_dev.data, b_dev.data);
	Info("initial residual: {}", sqrt(ArrayFunc::Dot(res_dev.data, res_dev.data)));

	VCycleMultigrid<T> solver;
	solver.Init_Poisson(poisson, 2, 2);

	//solver.Apply(x_dev.data, b_dev.data);
	//poisson.Residual(res_dev.data, x_dev.data, b_dev.data);
	//Info("mg residual: {}", sqrt(ArrayFunc::Dot(res_dev.data, res_dev.data)));

	//Info("run mg again:");

	//for (int k = 0; k < 100; k++) {
	//	solver.Apply(x_dev.data, b_dev.data);
	//	poisson.Residual(res_dev.data, x_dev.data, b_dev.data);
	//	Info("mg residual: {}", sqrt(ArrayFunc::Dot(res_dev.data, res_dev.data)));
	//}
}

template<class T, int d>
void Test_MGPCG(const Vector<int, d> counts) {
	Grid<d> grid(counts);
	PoissonMapping<T, d> poisson = Random_Poisson_Mapping<T, d>(grid);
	Field<T, d> b_host(grid);
	Random::Fill_Random_Array<T>(b_host.data, -5, 10);
	FieldDv<T, d> b_dev = b_host;
	FieldDv<T, d> x_dev(grid, 0);
	ConjugateGradient<T> MGPCG;
	VCycleMultigrid<T> precond;
	precond.Init_Poisson(poisson, 2, 2);
	MGPCG.Init(&poisson, &precond, false);
	//MGPCG.Init(&poisson, nullptr, true);
	int iters = 0;
	real res = 0;
	MGPCG.Solve(x_dev.data, b_dev.data, iters, res);
	Info("MGPCG solved {} iters with relative_error={}", iters, res);
}