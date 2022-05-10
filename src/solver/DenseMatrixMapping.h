//////////////////////////////////////////////////////////////////////////
// Dense matrix mapping
// Copyright (c) (2022-), Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "LinearMapping.h"
#include "AuxFunc.h"

namespace Meso {

	template<class T>
	class DenseMatrixMapping : LinearMapping<T> {
	public:
		int cols;//xdof
		int rows;//ydof
		ArrayDv<T> A;
		
		virtual int XDoF() const { return cols; }

		virtual int YDoF() const { return rows; }

		//template<int d>
		//void Init_PoissonLike(const Grid<d> grid, LinearMapping<T>& mapping) {
		//	cols = mapping.XDoF();
		//	rows = mapping.YDoF();
		//	Assert(cols == rows, "DenseMatrixMapping: cols={} mismatch rows={}", cols, rows);
		//	A.resize(cols * rows);
		//	temp_Ap.resize(cols);
		//	temp_p.resize(rows);
		//	ArrayFunc::Fill(A, (T)0);
		//	int row_nnz = (d == 2 ? 5 : 7);
		//	for (int flag = 0; flag < row_nnz; flag++) {//set all cells with color==flag to 1 and others to 0
		//		PoissonLikeMask<d> mask(flag);
		//		grid.Exec_Kernel(&Set_Cell_By_Color<T, d>, grid, mask, thrust::raw_pointer_cast(temp_p.data()));
		//		mapping.Apply(temp_Ap, temp_p);
		//		grid.Exec_Kernel(&Fill_Matrix_From_Result<T, d>, grid, mask, thrust::raw_pointer_cast(temp_Ap.data()), rows, thrust::raw_pointer_cast(A.data()));
		//	}
		//}

		//input p, get Ap
		virtual void Apply(ArrayDv<T>& Ap, const ArrayDv<T>& p) {
			Assert(false, "DenseMatrixMapping not implemented yet");
		}
	};
}