//////////////////////////////////////////////////////////////////////////
// Dense matrix mapping
// Copyright (c) (2022-), Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "LinearMapping.h"
#include "PoissonMapping.h"

namespace Meso {
	template<class T>
	class DenseMatrixMapping : LinearMapping<T> {
	public:
		int cols;//xdof
		int rows;//ydof
		ArrayDv<T> A;
		virtual int XDof() const { return cols; }

		virtual int YDof() const { return rows; }

		template<int d>
		void Init_PoissonLike(const Grid<d> grid, LinearMapping<T>& mapping) {
			cols = mapping.XDof();
			rows = mapping.YDof();
			Assert(cols == rows, "DenseMatrixMapping: cols={} mismatch rols={}", cols, rows);
			A.resize(cols * rows);
			int row_nnz = (d == 2 ? 5 : 7);
		}

		//input p, get Ap
		virtual void Apply(ArrayDv<T>& Ap, const ArrayDv<T>& p) {
			Assert(false, "DenseMatrixMapping not implemented yet");
		}
	};
}