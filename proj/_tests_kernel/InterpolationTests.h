//////////////////////////////////////////////////////////////////////////
// Test for interpolation
// Copyright (c) (2022-), Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "Interpolation.h"
#include "Random.h"
#include <limits>

using namespace Meso;

template<class T, int d>
void Test_Interpolation(const Vector<int, d> counts) {
	Typedef_VectorD(d);
	VectorD domain_min = MathFunc::V<d>(-0.9, -1.2, 3);
	Grid<d> grid(counts, 0.01, domain_min, CORNER);

	Field<T, d> my_field(grid);
	VectorD a = MathFunc::V<d>(0.2, 0.8, 0.3);
	T b = 7;
	auto f = [&](const VectorD pos) {return pos.dot(a) + b; };
	my_field.Calc_Nodes(
		[&](const VectorDi node) {
			return f(grid.Position(node));
		}
	);
	//VectorD domain_min = grid.Domain_Min();
	VectorD domain_max = domain_min + grid.dx * (counts - VectorDi::Ones()).template cast<real>();
	for (int i = 0; i < counts.prod() * 10; i++) {
		Vector<real, d> pos = Random::Uniform_In_Box(domain_min, domain_max);
		VectorDi coord; VectorD frac; grid.Get_Fraction(pos, coord, frac);
		Assert(grid.Valid(coord), "Test_Interpolation encountered an invalid position {} coord {} frac {} against counts {}", pos, coord, frac, grid.Counts());
		real from_intp = IntpLinear::Value(my_field, pos);
		//real from_intp = Interpolation<PointIntpLinear>::Value(grid, my_data.data(), coord, frac);
		if (fabs(from_intp - f(pos)) > MathFunc::Eps<T>()) {
			Info("coord: {}, frac: {}", coord, frac);
			Error("Test_Interpolation: counts={}, failed at {}", counts, pos);
			return;
		}
	}
	Pass("Test_Interpolation passed for counts={}", counts);
}