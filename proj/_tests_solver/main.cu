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
using namespace Meso;

int main(){
    Vector<int, 2> abc(1, 2);

    Test_Sparse_Matrix();
    Test_CG_Memory_Safe();

    Test_Poisson_Diagonal<float>(Vector2i(14, 53));
    Test_Poisson_Diagonal<double>(Vector2i(14, 53));
    Test_Poisson_Diagonal<float>(Vector3i(16, 44, 23));
    Test_Poisson_Diagonal<double>(Vector3i(16, 44, 23));

    //Test_Damped_Jacobian<float, 2>(10);//not a standard "yes-or-no" test
    //Test_Damped_Jacobian<double, 3>(10);//not a standard "yes-or-no" test

    Test_Coarsener2(Vector2i(14, 53));
    return 0;
}