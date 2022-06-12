//////////////////////////////////////////////////////////////////////////
// Smoother for Poisson Mapping
// Copyright (c) (2022-), Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "PoissonMapping.h"
#include "DenseMatrixMapping.h"
#include "SparseMatrixMapping.h"

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
	void PoissonLike_Diagonal(ArrayDv<T>& diag, MaskedPoissonMapping<T, d>& mapping) {
		const auto& grid = mapping.vol.grid;
		size_t n = mapping.XDoF();
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
		static constexpr int row_nnz = (d == 2 ? 5 : 7);
		int offset;
		PoissonLikeMask(const int _offset = 0) {
			offset = (_offset % row_nnz + row_nnz) % row_nnz;
		}
		__host__ __device__ int operator () (const VectorDi coord) const {
			if constexpr (d == 2) {
				return (coord[0] + coord[1] * 2 + offset) % 5;
			}
			else if constexpr (d == 3) {
				return (coord[0] + coord[1] * 2 + coord[2] * 3 + offset) % 7;
			}
			else Assert(false, "PoissonLikeMask not defined for d={}", d);
		}

		__host__ __device__ VectorDi Coord_Offset_To_Zero(const int mask_value) const {
			if constexpr (d == 2) {
				switch (mask_value) {
				case 0:return Vector2i(0, 0);
				case 1:return Vector2i(1, 0);
				case 4:return Vector2i(-1, 0);
				case 2:return Vector2i(0, 1);
				case 3:return Vector2i(0, -1);
				default:break;
				}
			}
			else if constexpr (d == 3) {
				switch (mask_value) {
				case 0:return Vector3i(0, 0, 0);
				case 1:return Vector3i(1, 0, 0);
				case 2:return Vector3i(0, 1, 0);
				case 3:return Vector3i(0, 0, 1);
				case 4:return Vector3i(0, 0, -1);
				case 5:return Vector3i(0, -1, 0);
				case 6:return Vector3i(-1, 0, 0);
				default:break;
				}
			}
			Assert(false, "PoissonLikeMask::Coord_Offset_To_Zero not defined for d={}", d);

		}
	};

	template<class T, int d>
	__global__ void Set_Cell_By_Color(const Grid<d> grid, const PoissonLikeMask<d> mask, T* cell_data) {
		Typedef_VectorD(d);
		VectorDi coord = GPUFunc::Thread_Coord<d>(blockIdx, threadIdx);
		int index = grid.Index(coord);
		if (mask(coord) == 0) cell_data[index] = 1;
		else cell_data[index] = 0;
	}
	//column-major
	template<class T, int d>
	__global__ void Fill_Dense_Matrix_From_Result(const Grid<d> grid, const PoissonLikeMask<d> mask, const T* Ap, const int ydof, T* mat) {
		Typedef_VectorD(d);
		VectorDi coord = GPUFunc::Thread_Coord<d>(blockIdx, threadIdx);
		//the cell that is switched to 1 for this time
		VectorDi on_coord = coord + mask.Coord_Offset_To_Zero(mask.row_nnz - mask(coord));
		if (grid.Valid(on_coord)) {
			int row_idx = grid.Index(coord), col_idx = grid.Index(on_coord);
			mat[row_idx + ydof * col_idx] = Ap[row_idx];
		}
	}
	//column-major
	template<class T, int d>
	void Dense_Matrix_From_PoissonLike(int& cols, int& rows, ArrayDv<T>& A, const Grid<d> grid, LinearMapping<T>& poisson_like, T diag_add_epsilon = 0) {
		ArrayDv<T> temp_p, temp_Ap;
		cols = poisson_like.XDoF();
		rows = poisson_like.YDoF();
		Assert(cols == rows, "Dense_Matrix_From_Poisson_Like: cols={} mismatch rows={}", cols, rows);
		A.resize(cols * rows);
		temp_Ap.resize(cols);
		temp_p.resize(rows);
		ArrayFunc::Fill(A, (T)0);
		int row_nnz = (d == 2 ? 5 : 7);
		for (int flag = 0; flag < row_nnz; flag++) {//set all cells with color==flag to 1 and others to 0
			PoissonLikeMask<d> mask(flag);
			grid.Exec_Kernel(&Set_Cell_By_Color<T, d>, grid, mask, ArrayFunc::Data<T, DEVICE>(temp_p));
			poisson_like.Apply(temp_Ap, temp_p);
			grid.Exec_Kernel(&Fill_Dense_Matrix_From_Result<T, d>, grid, mask, ArrayFunc::Data<T, DEVICE>(temp_Ap), rows, ArrayFunc::Data<T, DEVICE>(A));
		}
		for (int i = 0; i < rows; i++) {
			A[i * rows + i] += diag_add_epsilon;
		}
	}
	//column-major
	template<class T, int d>
	void DenseMatrixMapping_From_PoissonLike(DenseMatrixMapping<T>& dense_mapping, const Grid<d> grid, LinearMapping<T>& poisson_like, T diag_add_epsilon = 0) {
		Dense_Matrix_From_PoissonLike(dense_mapping.cols, dense_mapping.rows, dense_mapping.A, grid, poisson_like, diag_add_epsilon);
	}

	//Will add epsilon*I to the system
	//The reason of that feature is that a Poisson system may have a eigen value 0, add epsilon*I will make it positive definite
	template<class T, int d>
	SparseMatrix<T> SparseMatrix_From_PoissonLike(const Grid<d> grid, LinearMapping<T>& poisson_like, T diag_add_epsilon = 0) {
		int cols, rows;
		ArrayDv<T> A_dev;
		//column-major
		Dense_Matrix_From_PoissonLike(cols, rows, A_dev, grid, poisson_like, diag_add_epsilon);
		Array<T> A_host = A_dev;
		std::vector<Eigen::Triplet<T, int>> elements;
		//Info("rows {} cols {}", rows, cols);
		//fmt::print("[");
		for (int i = 0; i < rows; i++) {
			//fmt::print("[");
			for (int j = 0; j < cols; j++) {
				int idx = j * rows + i;
				T a = A_host[idx];
				if (a != 0) {
					elements.push_back(Eigen::Triplet<T, int>(i, j, a));
				}
				//fmt::print("{},", a);
			}
			//fmt::print("],\n");
		}
		//fmt::print("]");
		Eigen::SparseMatrix<T, Eigen::RowMajor, int> sparse_host;
		sparse_host.resize(rows, cols);
		sparse_host.setFromTriplets(elements.begin(), elements.end());
		sparse_host.makeCompressed();
		return sparse_host;
	}
}