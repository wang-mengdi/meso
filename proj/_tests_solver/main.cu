#include <iostream>
#include <cstdio>

#include "Common.h"
#include "Grid.h"
#include "Field.h"
#include <fmt/ranges.h>
#include "SparseTests.h"
#include "PoissonMapping.h"

int main(){
    /*PoissonMapping<real, 2> mapping;
    Array<int, DEVICE> a;
    Array<int, DEVICE> b;
    Array<int, DEVICE> c;
    auto plus = [=] __device__ (int i, int j)->int {return i + j; };
    ArrayFunc::Binary_Transform(a, b, plus, c);*/

    Test_Sparse_Matrix();
    

    return 0;
}