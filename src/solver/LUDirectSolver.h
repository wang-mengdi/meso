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

		virtual int XDof() const { return dof; }

		virtual int YDof() const { return dof; } 

		LUDenseSolver() {}
		LUDenseSolver(const DenseMatrixMapping<T>& dense_mapping) { Init(dense_mapping); }

		void Init(const DenseMatrixMapping<T>& dense_mapping) {
			dof = dense_mapping.XDof();
			Assert(dense_mapping.YDof() == dof, "LUDenseSolver::Init(): must input a square matrix");

			A = dense_mapping.A;//deep copy here
			piv.resize(dof);

			info.resize(1);

			T* A_ptr = thrust::raw_pointer_cast(A.data());
			int buffer_size;
			cusolverDnCreate(&solve_handle);
			checkCudaErrors(cudaGetLastError());

			if constexpr (std::is_same<T, double>::value)cusolverDnDgetrf_bufferSize(solve_handle, dof, dof, A_ptr, dof, &buffer_size);
			else cusolverDnSgetrf_bufferSize(solve_handle, dof, dof, A_ptr, dof, &buffer_size);
			buffer.resize(buffer_size);
			cudaDeviceSynchronize();
			checkCudaErrors(cudaGetLastError());

			T* buffer_ptr = thrust::raw_pointer_cast(buffer.data());
			int* piv_ptr = thrust::raw_pointer_cast(piv.data());
			ArrayDv<int> info(1);
			int* info_ptr = thrust::raw_pointer_cast(info.data());
			if constexpr (std::is_same<T, double>::value) cusolverDnDgetrf(solve_handle, dof, dof, A_ptr, dof, buffer_ptr, piv_ptr, info_ptr);
			else cusolverDnSgetrf(solve_handle, dof, dof, A_ptr, dof, buffer_ptr, piv_ptr, info_ptr);
			cudaDeviceSynchronize();
			checkCudaErrors(cudaGetLastError());
		}

		//input b, get x
		virtual void Apply(ArrayDv<T>& x, const ArrayDv<T>& b) {
			ArrayFunc::Copy(x, b);

			T* A_ptr = ArrayFunc::Data<T, DEVICE>(A);
			int* piv_ptr = ArrayFunc::Data<int, DEVICE>(piv);
			T* x_ptr = ArrayFunc::Data<T, DEVICE>(x);

			Info("A size: {}, piv size: {}, x size: {}", A.size(), piv.size(), x.size());

			int* info_ptr = ArrayFunc::Data<int, DEVICE>(info);
			if constexpr (std::is_same<T, double>::value) cusolverDnDgetrs(solve_handle, CUBLAS_OP_N, dof, 1, A_ptr, dof, piv_ptr, x_ptr, dof, info_ptr);
			else cusolverDnSgetrs(solve_handle, CUBLAS_OP_N, dof, 1, A_ptr, dof, piv_ptr, x_ptr, dof, info_ptr);
			cudaDeviceSynchronize();
			checkCudaErrors(cudaGetLastError());
		}
	};
}