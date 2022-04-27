//////////////////////////////////////////////////////////////////////////
// Grid Colored Gauss-Seidel Smoother
// Copyright (c) (2022-), Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "LinearMapping.h"
#include "PoissonFunc.h"
using namespace thrust::placeholders;

namespace Meso {
	template<int d>
	class GridGSMask {
		Typedef_VectorD(d);
	public:
		Grid<d> grid;
		GridGSMask() {};
		GridGSMask(const Grid<d> _grid) {
			grid = _grid;
		}
		__host__ __device__ int operator () (const int idx)const {
			VectorDi coord = grid.Coord(idx);
			int col = (coord[0] & 1);
			col |= ((coord[1] & 1) << 1);
			if (d == 3) col |= ((coord[2] & 1) << 2);
		}
	};

	template<class T, int d>
	class GridGSSmoother : public LinearMapping<T> {
	public:
		static constexpr int color_num = (d == 2 ? 4 : 8);
		LinearMapping<T>* mapping;
		int dof;
		int iter_num;
		ArrayDv<T> diag;
		ArrayDv<T> x_temp;
		GridGSMask<d> mask;
		GridGSSmoother() {}
		GridGSSmoother(MaskedPoissonMapping<T, d>& _mapping, const int _iter_num) {
			Init_Poisson(_mapping, _iter_num);
		}
		void Init_Poisson(MaskedPoissonMapping<T, d>& _mapping, const int _iter_num) {
			mapping = &_mapping;
			iter_num = _iter_num;
			dof = mapping->XDof();
			Poisson_Diagonal(diag, _mapping);
			x_temp.resize(dof);
			mask = GridGSMask<d>(_mapping.Grid());
		}
		virtual int XDof()const { return dof; }
		virtual int YDof()const { return dof; }
		virtual void Apply(ArrayDv<T>& x, const ArrayDv<T>& b) {
			Memory_Check(x, b, "GridGSSmoother::Apply error: not enough memory space");
			ArrayFunc::Fill(x, (T)0);
			if (iter_num == 0) return;
			thrust::counting_iterator<int> idxbegin(0);
			thrust::counting_iterator<int> idxend = idxbegin + dof;
			for (int iter = 0; iter < iter_num; iter++) {
				for (int c = 0; c < color_num; c++) {
					mapping->Residual(x_temp, x, b);
					ArrayFunc::Divide(x_temp, diag);
					thrust::transform_if(
						x.begin(),//first1
						x.end(),//last1
						x_temp.begin(),//first2
						idxbegin,//stencil
						x.begin(),//result
						_1 + _2,//binary op
						[=]__device__(const int idx)->bool { return (mask(idx) == c); }
					);
				}
			}
		}
	};
}