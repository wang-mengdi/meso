//////////////////////////////////////////////////////////////////////////
// Sparse matrix mapping
// Copyright (c) (2018-), Bo Zhu, Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "Common.h"
#include "AuxFunc.h"
#include "LinearMapping.h"
#include <cusparse.h>
//Note: this header should not include AuxFuncCuda.h

namespace Meso {

	template<class T>
	cusparseDnVecDescr_t Create_DnVecDescr_t(T* x, int size) {//x must be on device
		cusparseDnVecDescr_t vec_t = nullptr;
		cusparseStatus_t stat = cusparseCreateDnVec(&vec_t, size, (void*)x, GPUFunc::Cuda_Real_Type<T>());
		return vec_t;
	}


	////The sparse matrix is stored in CRS format
	template<class T, DataHolder side> class SparseMatrixMapping : public LinearMapping<T>
	{
	public:
		bool realloc_on_shrink = true;
		////crs matrix
		int m = 0;					////rows
		int n = 0;					////cols
		int nnz = 0;					////nonzeros
		Array<int, side> ptr, col;
		Array<T, side> val;
		cusparseHandle_t cusparseHandle = nullptr;

		SparseMatrixMapping(const Eigen::SparseMatrix<T, Eigen::RowMajor, int>& A) {
			cusparseCreate(&cusparseHandle);
			m = A.rows(); n = A.cols(); nnz = A.nonZeros();

			ptr.resize(m + 1);
			col.resize(nnz);
			val.resize(nnz);

			int* A_ptr = const_cast<int*>(A.outerIndexPtr());
			int* A_col = const_cast<int*>(A.innerIndexPtr());
			T* A_val = const_cast<T*>(A.valuePtr());

			thrust::copy(A_ptr, A_ptr + m + 1, ptr.begin());
			thrust::copy(A_col, A_col + nnz, col.begin());
			thrust::copy(A_val, A_val + nnz, val.begin());
		}
		~SparseMatrixMapping() {
			if (cusparseHandle) cusparseDestroy(cusparseHandle);
		}
		////Eigen SparseMatrix interfaces
		//int rows() const { return m; }
		//int cols() const { return n; }
		//int outerSize() const { return m; }
		//int nonZeros() const { return nnz; }
		constexpr T* valuePtr() { return ArrayFunc::Data<T, side>(val); }
		constexpr const T* valuePtr() const { return ArrayFunc::Data<T, side>(val); }
		constexpr int* outerIndexPtr() { return ArrayFunc::Data<int, side>(ptr); }
		constexpr const int* outIndexPtr() const { return ArrayFunc::Data<int, side>(ptr); }
		constexpr int* innerIndexPtr() { return ArrayFunc::Data<int, side>(col); }
		constexpr const int* innerIndexPtr() const { return ArrayFunc::Data<int, side>(col); }
		//similar to SparseMatrix::resize(). May realloc ptr, but not others. With default option, do not change data side.
		//void resize(int _m, int _n, enum DataHolder new_side = DataHolder::UNKNOWN);
		//void resizeNonZeros(int _nnz);//resize() must be called before, so we assume data side already set here

		////Cuda SparseMatrix interfaces
		cusparseSpMatDescr_t Get_SpMatCescr_t(void) {
			if constexpr (side == DataHolder::DEVICE) {
				cusparseSpMatDescr_t matA = nullptr;
				cusparseStatus_t status = cusparseCreateCsr(&matA, m, n, nnz,
					thrust::raw_pointer_cast(ptr.data()), thrust::raw_pointer_cast(col.data()), thrust::raw_pointer_cast(val.data()),
					CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
					CUSPARSE_INDEX_BASE_ZERO, GPUFunc::Cuda_Real_Type<T>());
				return matA;
			}
			else {
				return nullptr;
			}
		}

		virtual int XDoF() const { return n; }//number of cols

		virtual int YDoF() const { return m; }//number of rows

		//input p, get Ap
		virtual void Apply(ArrayDv<T>& Ap, const ArrayDv<T>& p) {
			if constexpr (side == DataHolder::DEVICE) {
				Ap.resize(YDoF());
				thrust::fill(Ap.begin(), Ap.end(), (T)0);

				T one = 1;
				T zero = 0;
				T* alpha = &one, * beta = &zero;
				auto data_type = GPUFunc::Cuda_Real_Type<T>();

				//A,x,y must be on device
				//This function will internally allocate and erase temporary space dBuffer
				cusparseSpMatDescr_t A_desc = Get_SpMatCescr_t();
				cusparseDnVecDescr_t x_desc = Create_DnVecDescr_t(thrust::raw_pointer_cast(p.data()), n);
				cusparseDnVecDescr_t y_desc = Create_DnVecDescr_t(thrust::raw_pointer_cast(Ap.data()), m);
				size_t buffersize; void* dBuffer = nullptr;
				cusparseSpMV_bufferSize(
					cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE,
					alpha, A_desc, x_desc, beta, y_desc, data_type,
					CUSPARSE_MV_ALG_DEFAULT, &buffersize);
				cudaMalloc(&dBuffer, buffersize);
				cusparseSpMV(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE,
					alpha, A_desc, x_desc, beta, y_desc, data_type,
					CUSPARSE_MV_ALG_DEFAULT, dBuffer);
				cudaFree(dBuffer);
			}
			else {
				Assert(false, "[LMSolver::SparseMatrixCPX::applyMapping] can't apply mapping on host");
			}
		}
	};

	template<class T>
	class SparseDiagonalPreconditioner : public LinearMapping<T> {
	public:
		ArrayDv<T> diag_inv;
		int cols = 0, rows = 0;
		SparseDiagonalPreconditioner() {}
		SparseDiagonalPreconditioner(const SparseMatrixMapping<T, DEVICE>& A) { Init(A); }
		//TODO: support multi time Init() here
		void Init(const SparseMatrixMapping<T, DEVICE>& A) {
			cols = A.XDoF();
			rows = A.YDoF();
			Assert(cols == rows, "DiagonalPreconditioner::Init column number doesn't equal to row number");
			Assert(cols > 0, "DiagonalPreconditioner::Init can't solve an empty system");

			const int* ptr_dev = A.outIndexPtr();
			const int* col_dev = A.innerIndexPtr();
			const T* val_dev = A.valuePtr();
			//note: if you write something like A_dev.mat_dev.ptr in the lambda function, it will fail
			//seems __device__ lambda can't properly capture object member
			auto f = [ptr_dev, col_dev, val_dev]__device__(const int i)->T {
				T inv_i;
				//row of i
				int row_start = ptr_dev[i], row_end = ptr_dev[i + 1];
				bool found = false;
				for (int it = row_start; it < row_end; it++) {
					int col_num = col_dev[it];
					if (col_num == i) {
						T val = val_dev[it];
						if (val != 0) {
							found = true;
							inv_i = 1.0 / val;
						}
						break;
					}
					else if (col_num > i) break;
				}
				if (!found) {
					inv_i = 1.0;
				}
				return inv_i;
			};

			thrust::counting_iterator<int> idx_begin(0);
			thrust::counting_iterator<int> idx_end = idx_begin + rows;
			diag_inv.resize(cols);
			thrust::transform(idx_begin, idx_end, diag_inv.begin(), f);
			//cwise_mapping_with_idx_wrapper(diag_inv_dev, f, row_num);
			checkCudaErrors(cudaGetLastError());
		}

		virtual int XDoF() const { return cols; }
		virtual int YDoF() const { return rows; }

		//input p, get Ap
		virtual void Apply(ArrayDv<T>& Ap, const ArrayDv<T>& p) {
			ArrayFunc::Binary_Transform(p, diag_inv, thrust::multiplies<T>(), Ap);
		}
	};


}