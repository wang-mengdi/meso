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
	template<class T, int d>
	__global__ void Fill_Matrix_From_Result(const Grid<d> grid, const PoissonLikeMask<d> mask, const T* Ap, const int ydof, T* mat) {
		Typedef_VectorD(d);
		VectorDi coord = VectorFunc::Vi<d>(
			blockIdx.x * grid.block_size + threadIdx.x,
			blockIdx.y * grid.block_size + threadIdx.y,
			blockIdx.z * grid.block_size + threadIdx.z
			);
		//the cell that is switched to 1 for this time
		VectorDi on_coord = coord + mask.Coord_Offset_To_Zero(mask.row_nnz - mask(coord));
		if (grid.Valid(on_coord)) {
			int row_idx = grid.Index(coord), col_idx = grid.Index(on_coord);
			mat[row_idx + ydof * col_idx] = Ap[row_idx];
		}
	}

	template<class T>
	class DenseMatrixMapping : LinearMapping<T> {
	public:
		int cols;//xdof
		int rows;//ydof
		ArrayDv<T> A;
		ArrayDv<T> temp_p;
		ArrayDv<T> temp_Ap;
		virtual int XDof() const { return cols; }

		virtual int YDof() const { return rows; }

		template<int d>
		void Init_PoissonLike(const Grid<d> grid, LinearMapping<T>& mapping) {
			cols = mapping.XDof();
			rows = mapping.YDof();
			Assert(cols == rows, "DenseMatrixMapping: cols={} mismatch rols={}", cols, rows);
			A.resize(cols * rows);
			temp_Ap.resize(cols);
			temp_p.resize(rows);
			ArrayFunc::Fill(A, (T)0);
			int row_nnz = (d == 2 ? 5 : 7);
			for (int flag = 0; flag < row_nnz; flag++) {//set all cells with color==flag to 1 and others to 0
				PoissonLikeMask<d> mask(flag);
				grid.Exec_Kernel(&Set_Cell_By_Color<T, d>, grid, mask, thrust::raw_pointer_cast(temp_p.data()));
				mapping.Apply(temp_Ap, temp_p);
				grid.Exec_Kernel(&Fill_Matrix_From_Result<T, d>, grid, mask, thrust::raw_pointer_cast(temp_Ap.data()), rows, thrust::raw_pointer_cast(A.data()));
			}
		}

		//input p, get Ap
		virtual void Apply(ArrayDv<T>& Ap, const ArrayDv<T>& p) {
			Assert(false, "DenseMatrixMapping not implemented yet");
		}
	};
}