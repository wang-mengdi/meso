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
	VectorD domain_min = VectorFunc::V<d>(-0.9, -1.2, 3);
	Grid<d, CORNER> grid(counts, 0.01, domain_min);
	//Grid<d, CORNER> grid(counts, 0.01);
	
	Array<T> my_data(grid.DoF());

	VectorD a = VectorFunc::V<d>(0.2, 0.8, 0.3);
	T b = 7;
	auto f = [&](const VectorD pos) {return pos.dot(a) + b; };
	grid.Exec_Nodes(
		[&](const VectorDi node) {
			VectorD pos = grid.Position(node);
			int index = grid.Index(node);
			my_data[index] = f(pos);
		}
	);
	//VectorD domain_min = grid.Domain_Min();
	VectorD domain_max = grid.Domain_Max();
	for (int i = 0; i < counts.prod() * 10; i++) {
		Vector<real, d> pos = Random::Uniform_In_Box(domain_min, domain_max);
		VectorDi coord; VectorD frac; grid.Get_Fraction(pos, coord, frac);
		Assert(grid.Valid(coord), "Test_Interpolation encountered an invalid position {}", pos);
		real from_intp = Interpolation::Linear_Intp(grid, my_data.data(), coord, frac);
		if (fabs(from_intp - f(pos)) > MathFunc::Eps<T>()) {
			Info("coord: {}, frac: {}", coord, frac);
			Error("Test_Interpolation: counts={}, failed at {}", counts, pos);
			return;
		}
	}
	Pass("Test_Interpolation passed for counts={}", counts);
}