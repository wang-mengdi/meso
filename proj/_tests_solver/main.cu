#include <iostream>
#include <cstdio>

#include "Common.h"
#include "Grid.h"
#include "Field.h"
#include <fmt/ranges.h>
#include "SparseTests.h"
#include "PoissonMapping.h"
#include "PoissonTests.h"
#include "OperatorTests.h"
#include "DenseSolverTests.h"
#include "RestrictorTests.h"
#include "ProlongatorTests.h"
#include "MultigridTests.h"
using namespace Meso;

int main(){
    Test_MGPCG<float>(Vector3i(128, 128, 128));
    //Test_Multigrid<double>(Vector3i(32, 64, 128));

    //Vector3i counts(32, 64, 128);
    //int L = 2;
    //Array<PoissonMapping<double, 3>> poissons;
    //poissons.resize(L + 1);
    //for (int i = 0; i <= 2; i++) {
    //    Grid<3> grid(counts);
    //    poissons[i].Init(grid);
    //    counts /= 2;
    //}

    return 0;

    Test_Sparse_Matrix();
    Test_CG_Memory_Safe();

    Test_Poisson_Diagonal<float>(Vector2i(14, 53));
    Test_Poisson_Diagonal<double>(Vector2i(14, 53));
    Test_Poisson_Diagonal<float>(Vector3i(16, 44, 23));
    Test_Poisson_Diagonal<double>(Vector3i(16, 44, 23));

    //Test_Damped_Jacobian<float, 2>(10);//not a standard "yes-or-no" test
    //Test_Damped_Jacobian<double, 3>(10);//not a standard "yes-or-no" test

    Test_Coarsener2(Vector2i(14, 53));
    Test_Coarsener3(Vector3i(16, 44, 23));

    Test_LU_Dense_Solver<float>(Vector2i(14, 13));
    Test_LU_Dense_Solver<double>(Vector2i(14, 13));
    Test_LU_Dense_Solver<float>(Vector3i(8, 7, 11));
    Test_LU_Dense_Solver<double>(Vector3i(8, 7, 11));

    Test_Restrictor<float>(Vector2i(14, 53));
    Test_Restrictor<double>(Vector2i(14, 53));
    Test_Restrictor<float>(Vector3i(15, 44, 23));
    Test_Restrictor<double>(Vector3i(15, 44, 23));

    Test_Prolongator<float>(Vector2i(14, 53));
    Test_Prolongator<double>(Vector2i(14, 53));
    Test_Prolongator<float>(Vector3i(15, 44, 23));
    Test_Prolongator<double>(Vector3i(15, 44, 23));

    Test_MGPCG<float>(Vector2i(64, 67));
    Test_MGPCG<double>(Vector2i(64, 67));
    Test_MGPCG<float>(Vector3i(34, 34, 66));
    Test_MGPCG<double>(Vector3i(34, 34, 66));
    return 0;
}