#include <iostream>
#include <cstdio>

#include "Common.h"
#include "Grid.h"
#include "Field.h"
#include <fmt/ranges.h>
#include "ConjugateGradient.h"
#include "SparseMatrixMapping.h"
#include "PoissonMapping.h"

using namespace Meso;

int main(){
    /*PoissonMapping<real, 2> mapping;
    Array<int, DEVICE> a;
    Array<int, DEVICE> b;
    Array<int, DEVICE> c;
    auto plus = [=] __device__ (int i, int j)->int {return i + j; };
    ArrayFunc::Binary_Transform(a, b, plus, c);*/


    // test the solve for a sparse matrix
    // rows and cols number of the matrx
    int rows = 10;
    int cols = 10;

    //create a diagonal dominant matrix
    std::vector<Eigen::Triplet<real>> tripletList;
    tripletList.reserve(3*rows); //three entries per row
    for (int i=0;i<rows;i++)
    {
        for (int j = std::max(i - 1,0); j <= std::min(i + 1,cols-1); j++) {
            if (i == j) { tripletList.push_back(Eigen::Triplet<real>(i, j, rand())); }
            else { tripletList.push_back(Eigen::Triplet<real>(i, j, (real)rand() / RAND_MAX)); }
        }
    }
    Eigen::SparseMatrix<real, Eigen::RowMajor, int> A(rows,cols);
    A.setFromTriplets(tripletList.begin(), tripletList.end());
    
    //create b through x to make sure a solution exists
    VectorXd x(cols);
    for (int i = 0; i < rows; i++) {
        x[i] = (real)rand() / RAND_MAX;
    }
    VectorXd b = A * x;
    std::cout << "A: \n" <<A.toDense() << std::endl;
    std::cout << "x:  " <<x.transpose() << std::endl;
    std::cout << "b:  " <<b.transpose() << std::endl;
    
    //Solve with Eigen
    Eigen::ConjugateGradient<SparseMatrix<real>, Eigen::Lower | Eigen::Upper> e_cg;
    e_cg.compute(A);
    x = e_cg.solve(b);
    std::cout << "Eigen CG solve iterations:     " << e_cg.iterations() << std::endl;
    std::cout << "Eigen CG solve estimated error: " << e_cg.error() << std::endl;
    std::cout << "Eigen CG solved x:" << x.transpose() << std::endl;
    
    //Solve with our CG solver with linear mapping
    SparseMatrixMapping<real, DEVICE> smm(A);
    ConjugateGradient<real> cg;
    cg.verbose = true;
    cg.Init(&smm, nullptr, e_cg.maxIterations(), e_cg.tolerance()); //use the same max iteration and tolerance as Eigen
    
    //Verify linear mapping first
    ArrayDv<real> x_d(cols), b_d(cols);
    for (int i = 0; i < cols; i++) { x_d[i] = x[i]; }
    smm.Apply(b_d,x_d);
    
    VectorXd x_cg(cols);
    for (int i = 0; i < cols; i++) { x_cg[i] = x_d[i]; }
    
    if (x_cg.isApprox(x)) {
        Pass("Correct SparseMatrixMapping!");
    }
    else {
        Error("Linear mapped Ap not equal to b, linear mapped Ap:");
        for (int i = 0; i < cols; i++) { std::cout << b_d[i] << " "; b_d[i] = b[i]; }
        std::cout << std::endl;
    }
    
    //Verify our CG solver
    int iters = 0;
    real relative_error = 0;
    cg.Solve(x_d, b_d, iters, relative_error);
    for (int i = 0; i < cols; i++) { x_cg[i] = x_d[i]; }

    if (x_cg.isApprox(x)) {
        Pass("Correct Result!");
    }
    else {
        Error("Incorrect Result!");
        std::cout <<"Our x:"<<x_cg.transpose() << std::endl;
    }

    return 0;
}