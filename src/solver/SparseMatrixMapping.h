//////////////////////////////////////////////////////////////////////////
// Sparse CUDA
// Copyright (c) (2018-), Bo Zhu, Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "Common.h"
#include "AuxFunc.h"
#include "LinearMapping.h"
#include <cusparse.h>
//Note: this header should not include AuxFuncCuda.h

template<class T>
cusparseDnVecDescr_t Create_DnVecDescr_t(T* x, int size) {//x must be on device
	cusparseDnVecDescr_t vec_t = nullptr;
	cusparseStatus_t stat = cusparseCreateDnVec(&vec_t, size, (void*)x, GPUFunc::Cuda_Real_Type<T>());
	return vec_t;
}

namespace LMSolver {
	////The sparse matrix is stored in CRS format
	template<class T, DataHolder side> class SparseMatrix: public LinearMapping<T>
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

		SparseMatrix(const Eigen::SparseMatrix<T, Eigen::RowMajor, int>& A) {
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
		~SparseMatrix() {
			if (cusparseHandle) cusparseDestroy(cusparseHandle);
		}
		////Eigen SparseMatrix interfaces
		//int rows() const { return m; }
		//int cols() const { return n; }
		//int outerSize() const { return m; }
		//int nonZeros() const { return nnz; }
		//T* valuePtr() { return val.data(); }
		//const T* valuePtr() const { return val.data(); }
		//int* outerIndexPtr() { return ptr.data(); }
		//const int* outIndexPtr() const { return ptr.data(); }
		//int* innerIndexPtr() { return col.data(); }
		//const int* innerIndexPtr() const { return col.data(); }
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

		virtual int xDoF() const { return n; }//number of cols

		virtual int yDoF() const { return m; }//number of rows

		//input p, get Ap
		virtual void applyMapping(ArrayDv<T>& Ap, const ArrayDv<T>& p) {
			if constexpr (side == DataHolder::DEVICE) {
				Ap.resize(yDoF());
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

}