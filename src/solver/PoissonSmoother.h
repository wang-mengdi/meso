//////////////////////////////////////////////////////////////////////////
// Smoother for Poisson Mapping
// Copyright (c) (2022-), Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "PoissonMapping.h"

namespace Meso {

	template<class T, int d, DataHolder side=DEVICE>
	void Poisson_Diagonal(const PoissonMapping<T, d>& mapping, Array<T, side>& diag) {
		const auto& grid = mapping.vol.grid;
		
	}

}