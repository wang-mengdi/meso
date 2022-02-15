#include <iostream>
#include <cstdio>

#include "Common.h"
#include "Grid.h"
#include "Field.h"
#include <fmt/ranges.h>
#include "ConjugateGradient.h"
#include "SparseMatrixMapping.h"


int main(){
    Matrix<real, 3> M;
    M << 1, 2, 3, 4, 5, 6, 7, 8, 9;
    auto M_sparse = M.sparseView();
    LMSolver::SparseMatrix<real, DataHolder::DEVICE> M_dev(M_sparse);
    return 0;
}