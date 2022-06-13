//////////////////////////////////////////////////////////////////////////
// Direct solver of sparse matrix
// Copyright (c) (2022-), Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include <cusolverSp.h>
#include "SparseMatrixMapping.h"

namespace Meso {
	template<class T>
	class CholeskySparseSolver :public LinearMapping<T> {
	public:
		using Base=LinearMapping<T>;
		SparseMatrixMapping<T, DEVICE> mapping;
		T tolerance = 0;
		cusparseMatDescr_t mat_descr;
		cusolverSpHandle_t solve_handle;
		//currently we don't have a default constructor of SparseMatrixMapping
		CholeskySparseSolver() {}
		//somehow you can't say _mapping is a const -- Mengdi
		CholeskySparseSolver(SparseMatrix<T>&& _mapping, T _tol = (T)0) : mapping(std::forward<SparseMatrix<T>>(_mapping)), tolerance(_tol) {
			cusolverSpCreate(&solve_handle);
			cusparseCreateMatDescr(&mat_descr);
			cusparseSetMatType(mat_descr, CUSPARSE_MATRIX_TYPE_GENERAL);
			cusparseSetMatIndexBase(mat_descr, CUSPARSE_INDEX_BASE_ZERO);
			Assert(XDoF() == YDoF(), "CholeskySparseSolver error: XDoF={} not equal to YDoF={}", XDoF(), YDoF());
			checkCudaErrors(cudaGetLastError());
		}
		~CholeskySparseSolver() {
			if (solve_handle) cusolverSpDestroy(solve_handle);
		}
		virtual int XDoF() const { return mapping.YDoF(); }
		virtual int YDoF() const { return mapping.XDoF(); }
		virtual void Apply(ArrayDv<T>& x, const ArrayDv<T>& b) {
			Base::Memory_Check(x, b, "CholeskySparseSolver::Apply failed: not enough memory space");
			int singularity = 0;
			if constexpr (std::is_same<T, double>::value) {
				cusolverSpDcsrlsvchol(
					solve_handle,
					mapping.m,
					mapping.nnz,
					mat_descr,
					mapping.valuePtr(),
					mapping.outIndexPtr(),
					mapping.innerIndexPtr(),
					ArrayFunc::Data(b),
					tolerance,
					/*reorder*/0,
					ArrayFunc::Data(x),
					&singularity
				);
			}
			else {
				cusolverSpScsrlsvchol(
					solve_handle,
					mapping.m,
					mapping.nnz,
					mat_descr,
					mapping.valuePtr(),
					mapping.outIndexPtr(),
					mapping.innerIndexPtr(),
					ArrayFunc::Data<T, DEVICE>(b),
					tolerance,
					/*reorder*/0,
					ArrayFunc::Data<T, DEVICE>(x),
					&singularity
				);
			}
			Assert(singularity == -1, "CholeskySparseSolver error: system not positive definite, singularity={}", singularity);
			checkCudaErrors(cudaGetLastError());
		}
	};
}
