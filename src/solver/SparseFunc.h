//////////////////////////////////////////////////////////////////////////
// Sparse solver
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include <Eigen/Sparse>
#include "Common.h"
#include "SparseMatrixMapping.h"

namespace Meso {
	namespace SparseFunc {
		////block matrix operations
		template<class T, int dim1, int dim2> void Add_Block(SparseMatrix<real>& K, const int K_i, const int K_j, const Matrix<T, dim1, dim2>& K_b, const int Kb_i = 0, const int Kb_j = 0)
		{
			for (int i = 0; i < dim1; i++)for (int j = 0; j < dim2; j++) { K.coeffRef(K_i * dim1 + i, K_j * dim2 + j) += K_b.coeff(Kb_i * dim1 + i, Kb_j * dim2 + j); }
		}

		template<int dim1, int dim2, class T_MAT> void Copy_Block(SparseMatrix<real>& K, const int K_i, const int K_j, const T_MAT& K_b, const int Kb_i = 0, const int Kb_j = 0)
		{
			for (int i = 0; i < dim1; i++)for (int j = 0; j < dim2; j++) { K.coeffRef(K_i * dim1 + i, K_j * dim2 + j) = K_b.coeff(Kb_i * dim1 + i, Kb_j * dim2 + j); }
		}

		template<int dim> void Set_Block(SparseMatrix<real>& K, const int K_i, const int K_j, const real value)
		{
			for (int i = 0; i < dim; i++)for (int j = 0; j < dim; j++) { K.coeffRef(K_i * dim + i, K_j * dim + j) = value; }
		}

		inline void Set_Value(SparseMatrix<real>& K, const real value)	/////set all nonzero values to be value but keep all indices
		{
			int nnz = (int)K.nonZeros();
#pragma omp parallel for
			for (int i = 0; i < nnz; i++) { K.valuePtr()[i] = value; }
		}
	};
};