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
#include "IOFunc.h"

namespace Meso {

	void Test_Coarsener2(const Vector2i counts);
	void Test_Coarsener3(const Vector3i counts);

	template<class T, int d>
	void Test_Multigrid(const Vector<int, d> counts) {
		Check_Cuda_Memory("test begin");

		Grid<d> grid(counts);
		MaskedPoissonMapping<T, d> poisson = Random_Poisson_Mapping<T, d>(grid, 10);

		Field<T, d> b_host(grid);
		Random::Fill_Random_Array<T>(b_host.Data(), -5, 10);
		FieldDv<T, d> b_dev = b_host;
		FieldDv<T, d> x_dev(grid, 0);

		Check_Cuda_Memory("allocated poisson");

		FieldDv<T, d> res_dev(grid);
		poisson.Residual(res_dev.Data(), x_dev.Data(), b_dev.Data());
		Info("initial residual: {}", sqrt(ArrayFunc::Dot(res_dev.Data(), res_dev.Data())));

		//VCycleMultigridIntp<T, d> solver;
		VCycleMultigrid<T, RestrictorSum<T, d>, ProlongatorSum<T, d>> solver;
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
		VCycleMultigridIntp<T, d> precond;
		//VCycleMultigridSum<T, d> precond;
		precond.Init_Poisson(poisson, 2, 2);
		//MGPCG.Init(&poisson, &precond, false);
		MGPCG.Init(&poisson, &precond, false, -1, 1e-6);
		//MGPCG.Init(&poisson, nullptr, true);
		Timer timer;
		auto [iters, relative_error] = MGPCG.Solve(x_dev.Data(), b_dev.Data());
		//Info("MGPCG solved {} iters with relative_error={}", iters, res);
		if (iters < 100) {
			Pass("MGPCG test passed in {}s for counts={}, with {} iters and relative_error={}", timer.Lap_Time(), counts, iters, relative_error);
		}
		else {
			Error("MGPCG test failed for counts={}, with {} iters and relative_error={}", counts, iters, relative_error);
		}
	}

	template<class T, int d>
	void Test_MGPCG_Dirichlet(const Vector<int, d> counts, bool output_x) {
		Typedef_VectorD(d);
		Grid<d> grid(counts, (real)100.0/counts[0]);

		//set A
		MaskedPoissonMapping<T, d> poisson(grid);
		Field<bool, d> fixed(grid);
		fixed.Iterate_Nodes([&](const VectorDi cell) {
			for (int axis = 0; axis < d; axis++)
				if (cell[axis] == 0 || cell[axis] == grid.Counts()[axis] - 1)
					return true;
			return false;
			});
		FaceField<T, d> vol(grid, 1);
		poisson.Init(fixed, vol);

		//set b
		Field<T, d> b_host(grid, T(2 * d));

		//move boundary value to rhs
		Field<T, d> delta(grid);
		delta.Calc_Nodes([&](const VectorDi cell) {
			for (int axis = 0; axis < d; axis++)
				if (cell[axis] == 0 || cell[axis] == grid.Counts()[axis] - 1) {
					VectorD pos = grid.Position(cell);
					real norm = pos.norm();
					return norm * norm;
				}
			return (real)0;
			});
		b_host += delta;
		
		//solve
		FieldDv<T, d> b_dev = b_host;
		FieldDv<T, d> x_dev(grid, 0);
		ConjugateGradient<T> MGPCG;
		VCycleMultigridIntp<T, d> precond;
		precond.Init_Poisson(poisson, 2, 2);
		MGPCG.Init(&poisson, &precond, false, -1, 1e-6);
		Timer timer;
		auto [iters, relative_error] = MGPCG.Solve(x_dev.Data(), b_dev.Data());
		Pass("MGPCG_Dirichlet test passed in {}s for counts={}, with {} iters and relative_error={}", timer.Lap_Time(), counts, iters, relative_error);
		
		//output x
		if (output_x) {
			Field<T, d> x_host = x_dev;
			std::string x_name = "x.vts";
			VTKFunc::Write_VTS(x_host, x_name);
		}
	}

	template<class T, int d>
	void Test_MGPCG_Dirichlet_Neumann(const Vector<int, d> counts, bool output_x) {
		Typedef_VectorD(d);
		Grid<d> grid(counts, (real)1.0 / counts[0]);

		//set A
		MaskedPoissonMapping<T, d> poisson(grid);
		Field<bool, d> fixed(grid);
		fixed.Iterate_Nodes([&](const VectorDi cell) {
			for (int axis = 0; axis < d; axis++)
				if ((cell[axis] == 0 || cell[axis] == grid.Counts()[axis] - 1) && (cell[1] != 0))
					return true;
			return false;
			});
		FaceField<T, d> vol(grid);
		vol.Calc_Faces([&](const int axis, const VectorDi face) {
			if (axis == 1 && face[1] == 0)
				return 0;
			return 1;
			});
		poisson.Init(fixed, vol);

		//set b
		Field<T, d> b_host(grid);
		b_host.Calc_Nodes([&](const VectorDi cell) {
			if (cell[1] == 0)
				return 1;
			return 0;
			});

		//solve
		FieldDv<T, d> b_dev = b_host;
		FieldDv<T, d> x_dev(grid, 0);
		ConjugateGradient<T> MGPCG;
		VCycleMultigridIntp<T, d> precond;
		precond.Init_Poisson(poisson, 2, 2);
		MGPCG.Init(&poisson, &precond, false, -1, 1e-6);
		Timer timer;
		auto [iters, relative_error] = MGPCG.Solve(x_dev.Data(), b_dev.Data());
		Pass("MGPCG_Dirichlet_Neumann test passed in {}s for counts={}, with {} iters and relative_error={}", timer.Lap_Time(), counts, iters, relative_error);

		//output x
		if (output_x) {
			Field<T, d> x_host = x_dev;
			std::string x_name = "x.vts";
			VTKFunc::Write_VTS(x_host, x_name);
		}
	}
}