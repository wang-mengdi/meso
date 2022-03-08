//////////////////////////////////////////////////////////////////////////
// Poisson Mapping Tests
// Copyright (c) (2022-), Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "PoissonSmoother.h"
#include "Random.h"

namespace Meso {

	template<class T, int d>
	void Test_Poisson_Diagonal(Vector<int, d> counts) {
		Typedef_VectorD(d);
		Grid<d> grid(counts);
		FaceField<T, d> vol(grid);
		Field<bool, d> fixed(grid);
		PoissonMapping<T, d> mapping;
		
		vol.Iterate_Faces(
			[&](const int axis, const VectorDi face) {
				return Random::Uniform(0, 1);
			}
		);
		fixed.Iterate_Cells(
			[&](const VectorDi cell) {
				return !(bool)Random::RandInt(0, 9);
			}
		);
		mapping.Init(Grid<d, CELL>(counts), vol, fixed);
	}

}