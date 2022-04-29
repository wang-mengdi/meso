#include <iostream>
#include <cstdio>

#include "Common.h"
#include "Grid.h"
#include "Field.h"
#include "GridTests.h"
#include "InterpolationTests.h"

using namespace Meso;

template<int d>
void Test_Grid(void) {
    Typedef_VectorD(d);
    Grid<2> grid(Vector2i(9, 7));
    Field<int, 2> F(grid);
    grid.Exec_Nodes(
        [&](const VectorDi& cell) {
            F(cell) = grid.Index(cell);
        }
    );
    Info("indicies: \n{}", F);
}

int main(){
    int n = 10;
    SparseMatrix<double> A = MatrixXd::Random(n, n).sparseView(0.5, 1);
    VectorXd b(n), x(n);
    Eigen::ConjugateGradient<SparseMatrix<double>, Eigen::Lower | Eigen::Upper, Eigen::IncompleteCholesky<double, Eigen::Lower | Eigen::Upper>> cg;
    //Eigen::ConjugateGradient<SparseMatrix<double>, Eigen::Lower, Eigen::IncompleteCholesky<double>> cg;
    cg.compute(A);
    x = cg.solve(b);
    x = cg.solve(b);
    return 0;

    //Eigen::ConjugateGradient<Eigen::SparseMatrix<double, Eigen::RowMajor, int>, Eigen::Lower | Eigen::Upper, Eigen::IncompleteCholesky<double>> cg;
    ////cg.setTolerance((real)1e-6);
    //Eigen::SparseMatrix<double, Eigen::RowMajor, int> A0;
    //cg.compute(A0);
    //return 0;

    //Test_Grid<2>();//a visual test

    Test_Grid_Index<float>(Vector2i(114, 514));
    Test_Grid_Index<double>(Vector2i(192, 168));
    Test_Grid_Index<float, 3>(Vector3i(1926, 8, 17));
    Test_Grid_Index<double>(Vector3i(62, 40, 21));

    Test_Face_Grid<2>(Vector2i(114, 514));
    Test_Face_Grid<3>(Vector3i(62, 40, 21));

    Test_Interpolation<float, 2>(Vector2i(114, 514));
    Test_Interpolation<double, 2>(Vector2i(192, 168));
    Test_Interpolation<float, 3>(Vector3i(16, 8, 17));
    Test_Interpolation<double, 3>(Vector3i(62, 40, 21));
    return 0;
}