#include <iostream>
#include <cstdio>

#include "Common.h"
#include "SparseTests.h"

using namespace Meso;

int main(){
    //Test_MGPCG<float>(Vector3i(512, 512, 512));
    //return 0;

    Test_Sparse_Matrix();
    Test_CG_Memory_Safe();

    //Test_Damped_Jacobian<float, 2>(10);//not a standard "yes-or-no" test
    //Test_Damped_Jacobian<double, 3>(10);//not a standard "yes-or-no" test

    return 0;
}