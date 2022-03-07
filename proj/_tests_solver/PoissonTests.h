//////////////////////////////////////////////////////////////////////////
// Poisson Mapping Tests
// Copyright (c) (2022-), Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#include "PoissonSmoother.h"

namespace Meso {

	template<class T, int d>
	void Test_Poisson_Diagonal(Vector<int, d> counts) {
		Grid<d> grid(counts);
		FaceField<T, d> vol(grid);
		PoissonMapping<T, d> mapping;
		
		//mapping.Init(Grid<d>(counts),)
	}

}