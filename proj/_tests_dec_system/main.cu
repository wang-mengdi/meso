#include <iostream>
#include <cstdio>

#include "Common.h"
#include "Grid.h"
#include "Field.h"
#include "SparseTests.h"
#include "PoissonMapping.h"
#include "PoissonTests.h"
#include "RestrictorTests.h"
#include "ProlongatorTests.h"
#include "MultigridTests.h"
#include "ExteriorCalculusTests.h"
#include "DenseSolverTests.h"
#include "SparseSolverTests.h"
#include "SmootherTests.h"
using namespace Meso;

int main(){
    //Test_GridGSSmoother<double>(Vector2i(16, 16));
    //Test_DampedJacobiSmoother<double>(Vector2i(16, 16));
    //Test_Multigrid<double>(Vector2i(14, 53));
    //Test_MGPCG<double>(Vector3i(128, 129, 128));
    //return 0;

    Test_LU_Dense_Solver<float>(Vector2i(14, 13));
    Test_LU_Dense_Solver<double>(Vector2i(14, 13));
    Test_LU_Dense_Solver<float>(Vector3i(8, 7, 11));
    Test_LU_Dense_Solver<double>(Vector3i(8, 7, 11));

    Test_Cholesky_Sparse_Solve<float>(Vector2i(14, 13));
    Test_Cholesky_Sparse_Solve<double>(Vector2i(14, 13));
    Test_Cholesky_Sparse_Solve<float>(Vector3i(8, 7, 11));
    Test_Cholesky_Sparse_Solve<double>(Vector3i(8, 7, 11));

    Test_Poisson_Diagonal<float>(Vector2i(14, 53));
    Test_Poisson_Diagonal<double>(Vector2i(14, 53));
    Test_Poisson_Diagonal<float>(Vector3i(16, 44, 23));
    Test_Poisson_Diagonal<double>(Vector3i(16, 44, 23));

    Test_Exterior_Derivative_Cell<float>(Vector2i(140, 530));
    Test_Exterior_Derivative_Cell<double>(Vector2i(64, 32));
    Test_Exterior_Derivative_Cell<float>(Vector3i(32, 44, 64));
    Test_Exterior_Derivative_Cell<double>(Vector3i(16, 44, 23));

    Test_Exterior_Derivative_Face<float>(Vector2i(140, 530));
    Test_Exterior_Derivative_Face<double>(Vector2i(64, 32));
    Test_Exterior_Derivative_Face<float>(Vector3i(32, 44, 64));
    Test_Exterior_Derivative_Face<double>(Vector3i(16, 44, 23));

    Test_Coarsener2(Vector2i(14, 53));
    Test_Coarsener3(Vector3i(16, 44, 23));

    Test_Restrictor_Intp<float>(Vector2i(14, 53));
    Test_Restrictor_Intp<double>(Vector2i(14, 53));
    Test_Restrictor_Intp<float>(Vector3i(15, 44, 23));
    Test_Restrictor_Intp<double>(Vector3i(15, 44, 23));

    Test_Prolongator_Intp<float>(Vector2i(14, 53));
    Test_Prolongator_Intp<double>(Vector2i(14, 53));
    Test_Prolongator_Intp<float>(Vector3i(15, 44, 23));
    Test_Prolongator_Intp<double>(Vector3i(15, 44, 23));

    //Test_GridGSSmoother<double, 2>(Vector2i(64, 64));

    Test_MGPCG<float>(Vector2i(512, 513));
    Test_MGPCG<double>(Vector2i(512, 513));
    Test_MGPCG<float>(Vector3i(256, 256, 256));
    Test_MGPCG<double>(Vector3i(128, 129, 128));
    //Test_MGPCG<float>(Vector3i(64, 64, 64));
    //Test_MGPCG<double>(Vector3i(128, 129, 128));
    return 0;
}