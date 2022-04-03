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
	class ChessboardMask {
	public:
		Typedef_VectorD(d);
		Grid<d> grid;
		const int color;
		ChessboardMask(const Grid<d>& _grid, const int _color) :grid(_grid), color(_color) {}
		__host__ __device__ int operator () (const int idx) {
			VectorDi coord = grid.Coord(idx);
			//note that XNOR is NOT*XOR, and ^0x1 equals to NOT
			//so ^color^0x1 means "==color"
			if constexpr (d == 2) {
				return (coord[0] ^ coord[1] ^ color ^ 0x1) & 0x1;
			}
			else if constexpr (d == 3) {
				return (coord[0] ^ coord[1] ^ coord[2] ^ color ^ 0x1) & 0x1;
			}
			else {
				Assert(false, "Meso::ChessboardMask undefined for d=={}", d);
				return false;
			}
		}
	};

	template<class T, int d >
	void Poisson_Diagonal(ArrayDv<T>& diag, PoissonMapping<T, d>& mapping) {
		const auto& grid = mapping.vol.grid;
		size_t n = mapping.XDof();
		diag.resize(n);
		thrust::fill(diag.begin(), diag.end(), 0);
		ArrayDv<T> p_temp(n);
		ArrayDv<T> Ap_temp(n);
		thrust::counting_iterator<int> idxbegin(0);
		thrust::counting_iterator<int> idxend = idxbegin + n;

		////white mask
		ChessboardMask<d> white_mask(grid, 0);
		thrust::transform(idxbegin, idxend, p_temp.begin(), white_mask);
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

	//a mask to distinguish the dense elements in a poisson system
	template<int d>
	class PoissonLikeMask {
		Typedef_VectorD(d);
	public:
		__host__ __device__ int operator () (const VectorDi coord) {
			if constexpr (d == 2) {
				return (coord[0] + coord[1] * 2) % 5;
			}
			else if constexpr (d == 3) {
				return (coord[0] + coord[1] * 2 + coord[2] * 3) % 7;
			}
			else Assert(false, "PoissonLikeMask not defined for d={}", d);
		}
	};
}