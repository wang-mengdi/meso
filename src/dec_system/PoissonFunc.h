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
	template<class T>
	__global__ void PoissonLike_One_Over_Diagonal_Kernel2(const Grid<2> _grid, const unsigned char* _cell_type, 
		T** _vol, T* _one_over_diag)
	{
		Typedef_VectorD(2);
		// calculate index
		VectorDi coord = GPUFunc::Thread_Coord<2>(blockIdx, threadIdx);
		const int idx = threadIdx.x;
		const int idy = threadIdx.y;
		const int global_id = _grid.Index(coord);

		// define shared memory
		__shared__ T shared_vol_x[9][8];
		__shared__ T shared_vol_y[8][9];

		// load data
		unsigned char type = _cell_type[global_id];

		shared_vol_x[idx][idy] = _vol[0][_grid.Face_Index(0, coord)];
		if (idx == 7)
			shared_vol_x[8][idy] = _vol[0][_grid.Face_Index(0, coord + VectorDi::Unit(0))];

		shared_vol_y[idx][idy] = _vol[1][_grid.Face_Index(1, coord)];
		if (idy == 7)
			shared_vol_y[idx][8] = _vol[1][_grid.Face_Index(1, coord + VectorDi::Unit(1))];

		__syncthreads();

		T result = 0;
		if (type != 1 && type != 2)
		{
			result += shared_vol_x[idx][idy];
			result += shared_vol_x[idx + 1][idy];
			result += shared_vol_y[idx][idy];
			result += shared_vol_y[idx][idy + 1];
			result = 1.0 / result;
		}
		else
			result = 1.0;
		_one_over_diag[global_id] = result;
	}


	template<class T>
	__global__ void PoissonLike_One_Over_Diagonal_Kernel3(const Grid<3> _grid, const unsigned char* _cell_type,
		T** _vol, T* _one_over_diag)
	{
		Typedef_VectorD(3);
		// calculate index
		VectorDi coord = GPUFunc::Thread_Coord<3>(blockIdx, threadIdx);
		const int idx = threadIdx.x;
		const int idy = threadIdx.y;
		const int idz = threadIdx.z;
		const int global_id = _grid.Index(coord);

		// define shared memory
		__shared__ T shared_vol_x[5][4][4];
		__shared__ T shared_vol_y[4][5][4];
		__shared__ T shared_vol_z[4][4][5];

		// load data
		unsigned char type = _cell_type[global_id];

		shared_vol_x[idx][idy][idz] = _vol[0][_grid.Face_Index(0, coord)];
		if (idx == 3)
			shared_vol_x[4][idy][idz] = _vol[0][_grid.Face_Index(0, coord + VectorDi::Unit(0))];

		shared_vol_y[idx][idy][idz] = _vol[1][_grid.Face_Index(1, coord)];
		if (idy == 3)
			shared_vol_y[idx][4][idz] = _vol[1][_grid.Face_Index(1, coord + VectorDi::Unit(1))];

		shared_vol_z[idx][idy][idz] = _vol[2][_grid.Face_Index(2, coord)];
		if (idz == 3)
			shared_vol_z[idx][idy][4] = _vol[2][_grid.Face_Index(2, coord + VectorDi::Unit(2))];

		__syncthreads();

		T result = 0;
		if (type != 1 && type != 2)
		{
			result += shared_vol_x[idx][idy][idz];
			result += shared_vol_x[idx + 1][idy][idz];
			result += shared_vol_y[idx][idy][idz];
			result += shared_vol_y[idx][idy + 1][idz];
			result += shared_vol_z[idx][idy][idz];
			result += shared_vol_z[idx][idy][idz + 1];
			result = 1.0 / result;
		}
		else
			result = 1.0;
		_one_over_diag[global_id] = result;
	}

	template<class T, int d>
	void PoissonLike_One_Over_Diagonal(ArrayDv<T>& _one_over_diag, MaskedPoissonMapping<T, d>& _mapping)
	{
		Grid<d> grid = _mapping.Grid();
		if constexpr (d == 2)
			grid.Exec_Kernel(PoissonLike_One_Over_Diagonal_Kernel2<T>, grid, _mapping.cell_type.Data_Ptr(),
				ArrayFunc::Data(_mapping.vol.face_data_ptr), ArrayFunc::Data(_one_over_diag));
		else
			grid.Exec_Kernel(PoissonLike_One_Over_Diagonal_Kernel3<T>, grid, _mapping.cell_type.Data_Ptr(),
				ArrayFunc::Data(_mapping.vol.face_data_ptr), ArrayFunc::Data(_one_over_diag));
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
	__global__ void Set_Cell_By_Color(const GridIndexer<d> grid, const PoissonLikeMask<d> mask, T* cell_data) {
		Typedef_VectorD(d);
		VectorDi coord = GPUFunc::Thread_Coord<d>(blockIdx, threadIdx);
		int index = grid.Index(coord);
		if (mask(coord) == 0) cell_data[index] = 1;
		else cell_data[index] = 0;
	}
	//column-major
	template<class T, int d>
	__global__ void Fill_Dense_Matrix_From_Result(const GridIndexer<d> grid, const PoissonLikeMask<d> mask, const T* Ap, const int ydof, T* mat) {
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
	void Dense_Matrix_From_PoissonLike(int& cols, int& rows, ArrayDv<T>& A, const GridIndexer<d> grid, LinearMapping<T>& poisson_like, T diag_add_epsilon = 0) {
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
			grid.Exec_Kernel(&Set_Cell_By_Color<T, d>, grid, mask, ArrayFunc::Data(temp_p));
			poisson_like.Apply(temp_Ap, temp_p);
			grid.Exec_Kernel(&Fill_Dense_Matrix_From_Result<T, d>, grid, mask, ArrayFunc::Data(temp_Ap), rows, ArrayFunc::Data(A));
		}
		for (int i = 0; i < rows; i++) {
			A[i * rows + i] += diag_add_epsilon;
		}
	}
	//column-major
	template<class T, int d>
	void DenseMatrixMapping_From_PoissonLike(DenseMatrixMapping<T>& dense_mapping, const GridIndexer<d> grid, LinearMapping<T>& poisson_like, T diag_add_epsilon = 0) {
		Dense_Matrix_From_PoissonLike(dense_mapping.cols, dense_mapping.rows, dense_mapping.A, grid, poisson_like, diag_add_epsilon);
	}

	//Will add epsilon*I to the system
	//The reason of that feature is that a Poisson system may have a eigen value 0, add epsilon*I will make it positive definite
	template<class T, int d>
	SparseMatrix<T> SparseMatrix_From_PoissonLike(const GridIndexer<d> grid, LinearMapping<T>& poisson_like, T diag_add_epsilon = 0) {
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