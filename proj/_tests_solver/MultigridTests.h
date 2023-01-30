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
	/*
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
	*/
	template<class T, int d>
	void Test_MGPCG_Dirichlet(const Vector<int, d> counts, bool output_x) {
		Typedef_VectorD(d);
		real domain_size = 1.0;
		Grid<d> grid(counts, domain_size / counts[0], -VectorD::Ones() * domain_size * 0.5);
		//set A
		MaskedPoissonMapping<T, d> poisson(grid);
		Field<unsigned char, d> cell_type(grid);
		cell_type.Calc_Nodes([&](const VectorDi cell) {
			if (grid.Is_Boundary_Cell(cell))
				return 1;
			return 0;
			});
		FaceField<T, d> vol(grid, 1);
		poisson.Init(cell_type, vol);

		//set b
		Field<T, d> b_host(grid, T(-2 * d * grid.dx * grid.dx));
		//move boundary value to rhs
		grid.Iterate_Nodes([&](const VectorDi cell) {
			if(grid.Is_Boundary_Cell(cell)){
				VectorD pos = grid.Position(cell);
				real norm = pos.norm();
				real val = norm * norm;
				b_host(cell) = val;
				return;
			}
			});
		grid.Iterate_Nodes([&](const VectorDi cell) {
			if(grid.Is_Boundary_Cell(cell)) {
					real val = b_host(cell);
					for (int axis = 0; axis < d; axis++)
					{
						VectorDi left_cell = cell - VectorDi::Unit(axis);
						if (grid.Valid(left_cell) && !grid.Is_Boundary_Cell(left_cell))
								b_host(left_cell) += val;
						VectorDi right_cell = cell + VectorDi::Unit(axis);
						if (grid.Valid(right_cell) && !grid.Is_Boundary_Cell(right_cell))
								b_host(right_cell) += val;
					}
				}});
		
		//solve
		FieldDv<T, d> b_dev = b_host;
		FieldDv<T, d> x_dev(grid, 0);
		ConjugateGradient<T, d> MGPCG;
		VCycleMultigridIntp<T, d> precond;
		precond.Init_Poisson(poisson);
		MGPCG.Init(&poisson, &precond, false, -1, 1e-5);
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
	/*
	template<class T, int d>
	void Test_MGPCG_Dirichlet_Neumann(const Vector<int, d> counts, bool output_x) {
		Typedef_VectorD(d);
		Grid<d> grid(counts, (real)CommonConstants::pi / counts[0]);

		auto f = [&](const VectorDi& cell) {
			VectorD pos = grid.Position(cell);
			return exp(pos[0]) * sin(pos[1]);
		};

		//set A
		MaskedPoissonMapping<T, d> poisson(grid);
		Field<bool, d> fixed(grid);
		fixed.Calc_Nodes([&](const VectorDi cell) {
			for (int axis = 0; axis < d; axis++)
				if ((cell[axis] == 0 || cell[axis] == grid.Counts()[axis] - 1))
					return true;
			return false;
			});
		FaceField<T, d> vol(grid);
		vol.Calc_Faces([&](const int axis, const VectorDi face) {
			if (axis == 0 && face[0] <= 1)
				return 0;
			return 1;
			});
		poisson.Init(fixed, vol);

		//set b
		Field<T, d> b_host(grid);
		b_host.Calc_Nodes([&](const VectorDi cell)->T {
			if (fixed(cell)) {
				return f(cell);
			}
			else if (cell[0] == 1) {
				VectorDi nb_cell = cell; nb_cell[0] = 0;
				return -(f(cell) - f(nb_cell)) / (grid.dx * grid.dx);
			}
			else return 0;
			});

		//solve
		FieldDv<T, d> b_dev = b_host;
		FieldDv<T, d> x_dev(grid, 0);
		ConjugateGradient<T> MGPCG;
		VCycleMultigridIntp<T, d> precond;
		precond.Init_Poisson(poisson, 2, 2);
		MGPCG.Init(&poisson, &precond, true, -1);
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


	//template<class T, int d>
	//void Test_MGPCG_Dirichlet_Neumann(const Vector<int, d> counts, bool output_x) {
	//	Typedef_VectorD(d);
	//	Grid<d> grid(counts, (real)1.0 / counts[0]);

	//	//set A
	//	MaskedPoissonMapping<T, d> poisson(grid);
	//	Field<bool, d> fixed(grid);
	//	fixed.Calc_Nodes([&](const VectorDi cell) {
	//		for (int axis = 0; axis < d; axis++)
	//			if ((cell[axis] == 0 || cell[axis] == grid.Counts()[axis] - 1) && (cell[1] != 0))
	//				return true;
	//		return false;
	//		});
	//	FaceField<T, d> vol(grid);
	//	vol.Calc_Faces([&](const int axis, const VectorDi face) {
	//		if (axis == 1 && face[1] == 0)
	//			return 0;
	//		return 1;
	//		});
	//	poisson.Init(fixed, vol);

	//	//set b
	//	Field<T, d> b_host(grid);
	//	b_host.Calc_Nodes([&](const VectorDi cell) {
	//		if (cell[1] == 0)
	//			return 1;
	//		return 0;
	//		});

	//	//solve
	//	FieldDv<T, d> b_dev = b_host;
	//	FieldDv<T, d> x_dev(grid, 0);
	//	ConjugateGradient<T> MGPCG;
	//	VCycleMultigridIntp<T, d> precond;
	//	precond.Init_Poisson(poisson, 2, 2);
	//	MGPCG.Init(&poisson, &precond, false, -1, 1e-5);
	//	Timer timer;
	//	auto [iters, relative_error] = MGPCG.Solve(x_dev.Data(), b_dev.Data());
	//	Pass("MGPCG_Dirichlet_Neumann test passed in {}s for counts={}, with {} iters and relative_error={}", timer.Lap_Time(), counts, iters, relative_error);

	//	//output x
	//	if (output_x) {
	//		Field<T, d> x_host = x_dev;
	//		std::string x_name = "x.vts";
	//		VTKFunc::Write_VTS(x_host, x_name);
	//	}
	//}
	*/
	template<class T, int d>
	void Test_MGPCG_Neumann(const Vector<int, d> counts, bool output_x) {
		Typedef_VectorD(d);
		Grid<d> grid(counts, (real)1.0 / counts[0]);
		//set A
		MaskedPoissonMapping<T, d> poisson(grid);
		Field<unsigned char, d> cell_type(grid, 0);
		FaceField<T, d> vol(grid);
		vol.Calc_Faces([&](const int axis, const VectorDi face) {
			if (face[axis] == 0 || face[axis] == grid.Counts()[axis])
				return 0;
			return 1;
			});
		poisson.Init(cell_type, vol);

		//set b
		Field<T, d> b_host(grid);
		b_host.Calc_Nodes([&](const VectorDi cell) {
			T val = 0;
			for (int axis = 0; axis < d; axis++)
			{
				if (cell[axis] == 0)
					val += 1;
				if (cell[axis] == grid.Counts()[axis] - 1)
					val -= 1;
			}
			return val;
			});

		//solve
		FieldDv<T, d> b_dev = b_host;
		FieldDv<T, d> x_dev(grid, 0);
		ConjugateGradient<T, d> MGPCG;
		VCycleMultigridIntp<T, d> precond;
		precond.Init_Poisson(poisson, 2, 20, 20);
		MGPCG.Init(&poisson, &precond, false, -1, 1e-5, true);
		Timer timer;
		auto [iters, relative_error] = MGPCG.Solve(x_dev.Data(), b_dev.Data());
		Pass("MGPCG_Neumann test passed in {}s for counts={}, with {} iters and relative_error={}", timer.Lap_Time(), counts, iters, relative_error);

		//output x
		if (output_x) {
			Field<T, d> x_host = x_dev;
			std::string x_name = "x.vts";
			VTKFunc::Write_VTS(x_host, x_name);
		}
	}
}