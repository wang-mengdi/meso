//////////////////////////////////////////////////////////////////////////
// Test sparse matrix
// Copyright (c) (2022-), Fan Feng
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#include "SparseMatrixMapping.h"
using namespace Meso;

void Sparse_Diagonal_Dominant_Matrix(int cols, int rows, Eigen::SparseMatrix<real, Eigen::RowMajor, int>& mat);
void Test_Sparse_Matrix(void);
void Test_CG_Memory_Safe(void);