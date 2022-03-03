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
	__host__ __device__ int Chessboard_Mask(const Grid<d> &grid, const int idx, const int color) {
		Typedef_VectorD(d);
		VectorDi coord = grid.Coord(idx);
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
	void Poisson_Diagonal(PoissonMapping<T, d>& mapping, Array<T, side>& diag) {
		const auto& grid = mapping.vol.grid;
		thrust::fill(diag.begin(), diag.end(), 0);
		size_t n = diag.end() - diag.begin();
		ArrayDv<T> p_temp(n);
		ArrayDv<T> Ap_temp(n);
		thrust::counting_iterator<int> idxbegin(0);
		thrust::counting_iterator<int> idxend = idxbegin + n;
		////white mask
		thrust::transform(idxbegin, idxend, p_temp, std::bind(Chessboard_Mask<d>, grid, std::placeholders::_1, 0));
		mapping.Apply(Ap_temp, p_temp);
		//Ap*.=p, masking out black cells
		ArrayFunc::Multiply(Ap_temp, p_temp);
		ArrayFunc::Add(diag, Ap_temp);
		////black mask
		//change p_temp from white to black
		ArrayFunc::Unary_Transform(p_temp, 1 - thrust::placeholders::_1, p_temp);
		mapping.Apply(Ap_temp, p_temp);
		//Ap*.=p, masking out white cells
		ArrayFunc::Multiply(Ap_temp, p_temp);
		ArrayFunc::Add(diag, Ap_temp);
	}

}