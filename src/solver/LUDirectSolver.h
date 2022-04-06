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
		ArrayDv<T> buffer;
		cusolverDnHandle_t solve_handle;

		virtual int XDof() const { return dof; }

		virtual int YDof() const { return dof; } 

		void Init(const DenseMatrixMapping<T>& dense_mapping) {
			dof = dense_mapping.XDof();
			Assert(dense_mapping.YDof() == dof, "LUDenseSolver::Init(): must input a square matrix");

			A.resize(dof * dof);
			piv.resize(dof);

			T* A_ptr = thrust::raw_pointer_cast(A.data());
			int buffer_size;
			cusolverDnCreate(&solve_handle);
			cusolverDnDgetrf_bufferSize(
				solve_handle,
				dof,
				dof,
				A_ptr,
				dof,
				&buffer_size
			);
			buffer.resize(buffer_size);

			T* buffer_ptr = thrust::raw_pointer_cast(buffer.data());
			int* piv_ptr = thrust::raw_pointer_cast(piv.data());
			ArrayDv<int> info(1);
			int* info_ptr = thrust::raw_pointer_cast(info.data());
			cusolverDnDgetrf(
				solve_handle,
				dof,
				dof,
				A_ptr,
				dof,
				buffer_ptr,
				piv_ptr,
				info_ptr
			);
		}

		//input b, get x
		virtual void Apply(ArrayDv<T>& x, const ArrayDv<T>& b) {
			ArrayFunc::Copy(x, b);
			ArrayDv<int> info(1);
			cusolverDnDgetrs(
				solve_handle,
				CUBLAS_OP_N,
				dof,
				1,
				thrust::raw_pointer_cast(A.data()),
				dof,
				thrust::raw_pointer_cast(piv.data()),
				thrust::raw_pointer_cast(x.data()),
				dof,
				thrust::raw_pointer_cast(info.data())
			);
		}
	};
}