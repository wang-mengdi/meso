//////////////////////////////////////////////////////////////////////////
// Dense matrix mapping
// Copyright (c) (2022-), Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "LinearMapping.h"
#include "PoissonMapping.h"
#include "AuxFunc.h"
#include "PoissonFunc.h"

namespace Meso {

	template<class T, int d>
	__global__ void Set_Cell_By_Color(const Grid<d> grid, const PoissonLikeMask<d> mask, T* cell_data) {
		Typedef_VectorD(d);
		VectorDi coord = VectorFunc::Vi<d>(
			blockIdx.x * grid.block_size + threadIdx.x,
			blockIdx.y * grid.block_size + threadIdx.y,
			blockIdx.z * grid.block_size + threadIdx.z
			);
		int index = grid.Index(coord);
		if (mask(coord) == 0) cell_data[index] = 1;
		else cell_data[index] = 0;
	}

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
			ArrayFunc::Fill(A, (T)0);
			int row_nnz = (d == 2 ? 5 : 7);
			for (int flag = 0; flag < row_nnz; flag++) {//set all cells with color==flag to 1 and others to 0

			}
		}

		//input p, get Ap
		virtual void Apply(ArrayDv<T>& Ap, const ArrayDv<T>& p) {
			Assert(false, "DenseMatrixMapping not implemented yet");
		}
	};
}