//////////////////////////////////////////////////////////////////////////
// LU Direct solver for dense matrix
// Copyright (c) (2022-), Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "DenseMatrixMapping.h"
#include <cusolverDn.h>

namespace Meso {
	template<class T>
	class LUDenseSolver :public LinearMapping<T> {
	public:
		int dof;

		ArrayDv<T> A;
		ArrayDv<int> piv;
		ArrayDv<int> info;//seems you can't destroy it
		ArrayDv<T> buffer;
		cusolverDnHandle_t solve_handle;

		virtual int XDoF() const { return dof; }

		virtual int YDoF() const { return dof; } 

		LUDenseSolver() {}
		LUDenseSolver(const DenseMatrixMapping<T>& dense_mapping) { Init(dense_mapping); }

		~LUDenseSolver() {
			cusolverDnDestroy(solve_handle);
		}

		void Init(const DenseMatrixMapping<T>& dense_mapping) {
			dof = dense_mapping.XDoF();
			Assert(dense_mapping.YDoF() == dof, "LUDenseSolver::Init(): must input a square matrix");

			A = dense_mapping.A;//deep copy here
			piv.resize(dof);

			info.resize(100);

			T* A_ptr = thrust::raw_pointer_cast(A.data());
			int buffer_size;
			cusolverDnCreate(&solve_handle);

			cusolverStatus_t flg;
			if constexpr (std::is_same<T, double>::value)flg = cusolverDnDgetrf_bufferSize(solve_handle, dof, dof, A_ptr, dof, &buffer_size);
			else flg = cusolverDnSgetrf_bufferSize(solve_handle, dof, dof, A_ptr, dof, &buffer_size);
			Assert(flg == CUSOLVER_STATUS_SUCCESS, "LUDenseSolver buffersize failed {}", flg);
			buffer.resize(buffer_size);

			T* buffer_ptr = thrust::raw_pointer_cast(buffer.data());
			int* piv_ptr = thrust::raw_pointer_cast(piv.data());
			int* info_ptr = thrust::raw_pointer_cast(info.data());
			if constexpr (std::is_same<T, double>::value) flg = cusolverDnDgetrf(solve_handle, dof, dof, A_ptr, dof, buffer_ptr, piv_ptr, info_ptr);
			else flg = cusolverDnSgetrf(solve_handle, dof, dof, A_ptr, dof, buffer_ptr, piv_ptr, info_ptr);
			Assert(flg == CUSOLVER_STATUS_SUCCESS, "LUDenseSolver getrf failed {}", flg);
			checkCudaErrors(cudaGetLastError());

			//Info("matrix A:");
			//for (int i = 0; i < dof; i++) {
			//	for (int j = 0; j < dof; j++) {
			//		T ele = A[i + dof * j];
			//		fmt::print("{} ", ele);
			//	}
			//	fmt::print("\n");
			//}
		}

		//input b, get x
		virtual void Apply(ArrayDv<T>& x, const ArrayDv<T>& b) {
			Memory_Check(x, b, "LUDenseSolver::Apply failed: not enough memory space");

			Assert(ArrayFunc::Is_Finite<T, DEVICE>(b), "LUDirect solver failed: b={}", b);

			ArrayFunc::Copy(x, b);

			T* A_ptr = ArrayFunc::Data<T, DEVICE>(A);
			int* piv_ptr = ArrayFunc::Data<int, DEVICE>(piv);
			T* x_ptr = ArrayFunc::Data<T, DEVICE>(x);
			cusolverStatus_t flg;
			//int* info_ptr = ArrayFunc::Data<int, DEVICE>(info);
			int* info_ptr = thrust::raw_pointer_cast(info.data());

			if constexpr (std::is_same<T, double>::value) flg = cusolverDnDgetrs(solve_handle, CUBLAS_OP_N, dof, 1, A_ptr, dof, piv_ptr, x_ptr, dof, info_ptr);
			else flg = cusolverDnSgetrs(solve_handle, CUBLAS_OP_N, dof, 1, A_ptr, dof, piv_ptr, x_ptr, dof, info_ptr);
			Assert(flg == CUSOLVER_STATUS_SUCCESS, "LUDenseSolver solve failed {}", flg);
			checkCudaErrors(cudaGetLastError());

			Assert(ArrayFunc::Is_Finite<T, DEVICE>(x), "LUDirect solver failed: x={}", x);
			//Info("direct solve input:\n{}\noutput:\n{}", b, x);
		}
	};
}