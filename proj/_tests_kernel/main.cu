#include <iostream>
#include <cstdio>

#include "Common.h"
#include "Grid.h"
#include "Field.h"
#include "GridTests.h"
#include "InterpolationTests.h"
#include "MeshTests.h"
#include "ImplicitManifold.h"
#include "SparseTests.h"

using namespace Meso;


void Print_Grid_Index2(Vector2i counts) {
    Grid<2> grid(counts, 1.0 / 7, Vector2::Zero(), CORNER);
    Field<int, 2> F(grid);
    grid.Exec_Nodes(
        [&](const Vector2i& cell) {
            F(cell) = grid.Index(cell);
        }
    );
    Info("grid: {}", F.grid);
    Info("indicies: \n{}", F);
}

int main() {
    ////Part 1: test reservoir

    //a visual test, not included in general tests
    //Print_Grid_Index2(Vector2i(8, 8));
    //Print_Grid_Index2(Vector2i(64, 64));

    //Test for an unpadded grid, if coord-index is a one-to-one mapping
    //So here the dimensions must be unpadded
    Test_Grid_Index<float>(Vector2i(120, 520));
    Test_Grid_Index<double>(Vector2i(192, 168));
    Test_Grid_Index<float, 3>(Vector3i(1928, 8, 20));
    Test_Grid_Index<double>(Vector3i(64, 40, 24));

    Test_Face_Grid<2>(Vector2i(114, 514));
    Test_Face_Grid<3>(Vector3i(62, 40, 21));

    Test_Interpolation<float, 2>(Vector2i(114, 514));
    Test_Interpolation<double, 2>(Vector2i(192, 168));
    Test_Interpolation<float, 3>(Vector3i(16, 8, 17));
    Test_Interpolation<double, 3>(Vector3i(62, 40, 21));

    Test_Interpolation_Clamp<double, 2>(Vector2i(32, 32));
}