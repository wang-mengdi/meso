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

template<class T, int d>
void Test_Interpolation_Clamp(const Vector<int, d> counts) {
	Typedef_VectorD(d);
	VectorD domain_min = MathFunc::V<d>(0, 0, 0);
	Grid<d> grid(counts, (real)1/(counts[0]-1), domain_min, CORNER);

	Field<T, d> my_field(grid);
	auto f = [&](const VectorD pos) {return (pos-grid.Center()).norm(); };
	my_field.Calc_Nodes(
		[&](const VectorDi node) {
			return f(grid.Position(node));
		}
	);

	VectorD domain_max = domain_min + grid.dx * (counts - VectorDi::Ones()).template cast<real>();
	for (int axis=0;axis<d;axis++){
		for (int i = 0; i < counts.prod() ; i++) {
			VectorD mask = VectorD::Ones();
			mask[axis] = 0;
			VectorD border_min = domain_min - (real) 2 * grid.dx *VectorD::Unit(axis); //The border that is outside of the domain
			VectorD border_max = domain_min + (domain_max - domain_min).cwiseProduct(mask);
			VectorD pos = Random::Uniform_In_Box(border_min, border_max);
			VectorDi coord; VectorD frac; grid.Get_Fraction(pos, coord, frac);
			real from_intp = IntpLinearClamp::Value(my_field, pos);

			Vector<real, d> pos1 = pos;
			pos1[axis] += (domain_max[axis] - domain_min[axis] - (real)2*pos[axis]);
			VectorDi coord1; VectorD frac1; grid.Get_Fraction(pos1, coord1, frac1);
			real from_intp1 = IntpLinearClamp::Value(my_field, pos1);

			if (fabs(from_intp - from_intp1) > MathFunc::Eps<T>()) {
				Info("pos: {}, from_intp: {}, pos1: {}, from_intp1: {}", pos, from_intp, pos1, from_intp1);
				Error("Test_Interpolation_Clamp: not symmetric upon clamp on axis={}", axis);
				return;
			}
		}
	}

	Pass("Test_Interpolation_Clamp passed for counts={}", counts);
}