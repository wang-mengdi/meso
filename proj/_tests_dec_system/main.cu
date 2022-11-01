#include <iostream>
#include <cstdio>

#include "Common.h"
#include "Grid.h"
#include "Field.h"
#include "PoissonMapping.h"
#include "PoissonTests.h"
#include "RestrictorTests.h"
#include "ProlongatorTests.h"
#include "MultigridTests.h"
#include "ExteriorDerivativeTests.h"
#include "DenseSolverTests.h"
#include "SparseSolverTests.h"
#include "SmootherTests.h"
using namespace Meso;

int main(){
    Test_Exterior_Derivative_Cell<float, 2, ExteriorDerivativePadding0>(Vector2i(140, 530));
    Test_Exterior_Derivative_Cell<double, 2, ExteriorDerivativePadding0>(Vector2i(64, 37));
    Test_Exterior_Derivative_Cell<float, 3, ExteriorDerivativePadding0>(Vector3i(32, 47, 64));
    Test_Exterior_Derivative_Cell<double, 3, ExteriorDerivativePadding0>(Vector3i(16, 44, 23));

    Test_Exterior_Derivative_Face<float, 2, ExteriorDerivativePadding0>(Vector2i(140, 530));
    Test_Exterior_Derivative_Face<double, 2, ExteriorDerivativePadding0>(Vector2i(64, 37));
    Test_Exterior_Derivative_Face<float, 3, ExteriorDerivativePadding0>(Vector3i(32, 47, 64));
    Test_Exterior_Derivative_Face<double, 3, ExteriorDerivativePadding0>(Vector3i(16, 44, 23));

    //NOTE: these grid sizes must have zero padding
    Test_LU_Dense_Solver<float>(Vector2i(16, 16));
    Test_LU_Dense_Solver<double>(Vector2i(16, 16));
    Test_LU_Dense_Solver<float>(Vector3i(8, 8, 12));
    Test_LU_Dense_Solver<double>(Vector3i(8, 8, 12));

    Test_Cholesky_Sparse_Solve<float>(Vector2i(16, 16));
    Test_Cholesky_Sparse_Solve<double>(Vector2i(16, 16));
    Test_Cholesky_Sparse_Solve<float>(Vector3i(8, 8, 12));
    Test_Cholesky_Sparse_Solve<double>(Vector3i(8, 8, 12));

    Test_Poisson_Diagonal<float>(Vector2i(16, 56));
    Test_Poisson_Diagonal<double>(Vector2i(16, 56));
    Test_Poisson_Diagonal<float>(Vector3i(16, 44, 24));
    Test_Poisson_Diagonal<double>(Vector3i(16, 44, 24));

    Test_Coarsener2(Vector2i(16, 56));//14,53 padded
    Test_Coarsener3(Vector3i(16, 44, 24));///16,44,23 padded

    Test_Restrictor_Intp<float>(Vector2i(16, 56));//14,53 padded
    Test_Restrictor_Intp<double>(Vector2i(16, 56));
    Test_Restrictor_Intp<float>(Vector3i(16, 44, 24));//15,44,23 padded
    Test_Restrictor_Intp<double>(Vector3i(16, 44, 24));

    Test_Prolongator_Intp<float>(Vector2i(16, 56));//14,53 padded
    Test_Prolongator_Intp<double>(Vector2i(16, 56));
    Test_Prolongator_Intp<float>(Vector3i(16, 44, 24));//15,44,23 padded
    Test_Prolongator_Intp<double>(Vector3i(16, 44, 24));

    ////A visual test
    //Test_GridGSSmoother<double, 2>(Vector2i(64, 64));
    
    Test_MGPCG<float>(Vector2i(512, 520));//512,513 padded
    Test_MGPCG<double>(Vector2i(512, 520));
    Test_MGPCG<float>(Vector3i(256, 256, 256));
    Test_MGPCG<double>(Vector3i(128, 132, 128));//128,129,128 padded

    Test_MGPCG_Dirichlet<double>(Vector2i(256, 256), false);
    Test_MGPCG_Dirichlet<double>(Vector2i(512, 512), false);
    Test_MGPCG_Dirichlet<double>(Vector3i(128, 128, 128), false);
    Test_MGPCG_Dirichlet<double>(Vector3i(256, 256, 256), false);

    Test_MGPCG_Dirichlet_Neumann<double>(Vector2i(256, 256), false);
    Test_MGPCG_Dirichlet_Neumann<double>(Vector2i(512, 512), false);
    Test_MGPCG_Dirichlet_Neumann<double>(Vector3i(128, 128, 128), false);
    Test_MGPCG_Dirichlet_Neumann<double>(Vector3i(256, 256, 256), false);

    Test_MGPCG_Neumann<double>(Vector2i(256, 256), false);
    Test_MGPCG_Neumann<double>(Vector2i(512, 512), false);
    Test_MGPCG_Neumann<double>(Vector3i(256, 256, 256), false);
    return 0;
}