#include "SparseTests.h"
#include "ConjugateGradient.h"
#include "Random.h"

template<class T>
void Random_Sparse_Diagonal_Dominant_Matrix(int cols, int rows, Eigen::SparseMatrix<T, Eigen::RowMajor, int>& mat) {
    std::vector<Eigen::Triplet<T>> tripletList;
    tripletList.reserve(3 * rows); //three entries per row
    for (int i = 0; i < rows; i++)
    {
        for (int j = std::max(i - 1, 0); j <= std::min(i + 1, cols - 1); j++) {
            if (i == j) { tripletList.push_back(Eigen::Triplet<T>(i, j, Random::Random() * 100)); }
            else { tripletList.push_back(Eigen::Triplet<T>(i, j, Random::Random())); }
        }
    }
    mat.resize(rows, cols);
    mat.setFromTriplets(tripletList.begin(), tripletList.end());
}

void Test_Sparse_Matrix(void)
{
    // test the solve for a sparse matrix
    // rows and cols number of the matrx
    int rows = 10;
    int cols = 10;

    //create a diagonal dominant matrix
    Eigen::SparseMatrix<real, Eigen::RowMajor, int> A;
    Random_Sparse_Diagonal_Dominant_Matrix(rows, cols, A);

    //create b through x to make sure a solution exists
    VectorXd x(cols);
    for (int i = 0; i < rows; i++) {
        x[i] = Random::Random();
    }
    VectorXd b = A * x;
    //std::cout << "A: \n" << A.toDense() << std::endl;
    //std::cout << "x:  " << x.transpose() << std::endl;
    //std::cout << "b:  " << b.transpose() << std::endl;

    //Solve with Eigen
    Eigen::ConjugateGradient<SparseMatrix<real>, Eigen::Lower | Eigen::Upper, Eigen::IdentityPreconditioner> e_cg;
    e_cg.compute(A);
    x = e_cg.solve(b);
    //std::cout << "Eigen CG solve iterations:     " << e_cg.iterations() << std::endl;
    //std::cout << "Eigen CG solve estimated error: " << e_cg.error() << std::endl;
    //std::cout << "Eigen CG solved x:" << x.transpose() << std::endl;

    //Solve with our CG solver with linear mapping
    SparseMatrixMapping<real, DEVICE> smm(A);
    ConjugateGradient<real> cg;
    cg.verbose = true;
    cg.Init(&smm, nullptr, e_cg.maxIterations(), e_cg.tolerance()); //use the same max iteration and tolerance as Eigen

    //Verify linear mapping first
    ArrayDv<real> x_d(cols), b_d(cols);
    for (int i = 0; i < cols; i++) { x_d[i] = x[i]; }
    smm.Apply(b_d, x_d);

    VectorXd x_cg(cols);
    for (int i = 0; i < cols; i++) { x_cg[i] = x_d[i]; }

    Assert(x_cg.isApprox(x), "Test_Sparse_Matrix: sparse mapped Ap not equal to b");

    ArrayFunc::Copy(b_d, b);

    //Verify our CG solver
    int iters = 0;
    real relative_error = 0;
    cg.Solve(x_d, b_d, iters, relative_error);
    for (int i = 0; i < cols; i++) { x_cg[i] = x_d[i]; }

    if (x_cg.isApprox(x)) {
        Pass("Test_Sparse_Matrix passed");
    }
    else {
        Error("Incorrect Result!");
        std::cout << "Our x:" << x_cg.transpose() << std::endl;
    }
}

void Test_CG_Memory_Safe(void) {
    int rows = 10;
    int cols = 10;

    //create a diagonal dominant matrix
    Eigen::SparseMatrix<real, Eigen::RowMajor, int> A;
    Random_Sparse_Diagonal_Dominant_Matrix(rows, cols, A);

    SparseMatrixMapping<real, DEVICE> smm_A(A);
    ConjugateGradient<real> cg;
    cg.verbose = true;
    cg.Init(&smm_A, nullptr, 1000, 1e-6);

    Eigen::SparseMatrix<real, Eigen::RowMajor, int> B;
    Random_Sparse_Diagonal_Dominant_Matrix(rows, cols, A);
    SparseMatrixMapping<real, DEVICE> smm_B(B);
    cg.Init(&smm_B, nullptr, 1000, 1e-6);
    Pass("Passed initializing memory for multiple times!");
}