#include <iostream>
#include <cstdio>

#include "Common.h"
#include "Grid.h"
#include "Field.h"
#include "GridTests.h"
#include "InterpolationTests.h"

#include <fmt/ranges.h>

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
    //for (int i = 0; i < grid.counts[0]; i++) {
    //    for (int j = 0; j < grid.counts[1]; j++) {
    //        fmt::print("{} ", F(Vector2i(i, j)));
    //    }
    //    fmt::print("\n");
    //}
}

int main(){
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