//////////////////////////////////////////////////////////////////////////
// Smoother for Poisson Mapping
// Copyright (c) (2022-), Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "PoissonMapping.h"

namespace Meso {
	//color==0: white
	//color==1: black
	template<int d>
	__host__ __device__ int Chessboard_Mask(const Vector<int, d> coord, int color) {
		//note that XNOR is NOT*XOR, and ^0x1 equals to NOT
		//so ^color^0x1 means "==color"
		if constexpr (d == 2) {
			return (coord[0] ^ coord[1] ^ color ^ 0x1) & 0x1;
		}
		else if constexpr(d == 3) {
			return (coord[0] ^ coord[1] ^ coord[2] ^ color ^ 0x1) & 0x1;
		}
		else {
			Assert("Meso::Is_Chessboard_Mask undefined for d=={}", d);
			return false;
		}
	}

	template<class T, int d, DataHolder side=DEVICE>
	void Poisson_Diagonal(const PoissonMapping<T, d>& mapping, Array<T, side>& diag) {
		const auto& grid = mapping.vol.grid;
		thrust::fill(diag.begin(), diag.end(), 0);
	}

}