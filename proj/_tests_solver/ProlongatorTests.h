//////////////////////////////////////////////////////////////////////////
// Test prolongator
// Copyright (c) (2022-), Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "Prolongator.h"
#include "Restrictor.h"
#include "Field.h"
#include "Random.h"
using namespace Meso;

template<class T,int d>
void Test_Prolongator(const Vector<int, d> fine_counts) {
	//test that restrict(prolongate())==identity
	Typedef_VectorD(d);
	Grid<d> fine_grid(fine_counts, 0.5, VectorD::Zero());
	Grid<d> coarse_grid(fine_counts / 2, 1.0, VectorD::Zero());
	Field<T, d> coarse_host(coarse_grid, 0);
	coarse_host(VectorFunc::Vi<d>(5, 5, 5)) = 1;
	//Random::Fill_Random_Array<T>(coarse_host.data, -5.0, 10.0);
	Grid<d> coarse_pad_grid(coarse_grid.counts + VectorDi::Ones(), 1.0, -VectorD::Ones() * 1.0);
	Field<T, d> coarse_pad(coarse_pad_grid, 0);
	coarse_pad.Calc_Cells(
		[&](const VectorDi cell) ->T {
			VectorDi offset_cell = cell - VectorDi::Ones();
			if (coarse_grid.Valid(offset_cell)) return coarse_host(offset_cell);
			else return 0;
		}
	);

	FieldDv<T, d> coarse_dev = coarse_host;
	FieldDv<T, d> fine_dev(fine_grid);
	Prolongator<T, d> P; P.Init(fine_grid, coarse_grid);
	P.Apply(fine_dev.Data(), coarse_dev.Data());

	Field<T, d> fine_cpu(fine_grid);
	fine_cpu.Calc_Cells(
		[&](const VectorDi cell) {
			VectorD pos = fine_grid.Position(cell);
			return IntpLinear::Value(coarse_pad, pos);
			//return Interpolation::Linear_Intp(coarse_pad, pos);
		}
	);

	Field<T, d> result_host = fine_dev;
	//Info("fine_cpu: \n{}\n", fine_cpu);
	//Info("result_host: \n{}\n", result_host);
	if (ArrayFunc::IsApprox<T>(result_host.Data(), fine_cpu.Data())) {
		Pass("Test_Prolongator passed for fine_counts={}", fine_counts);
	}
	else {
		Error("Test_Prolongator failed for fine_counts={}", fine_counts);
	}
}
