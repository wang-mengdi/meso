#include <iostream>
#include <cstdio>

#include "Common.h"
#include "Grid.h"
#include "Field.h"
#include "GridTests.h"
#include "InterpolationTests.h"
#include "MeshTests.h"

#include "SparseTests.h"

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

int main() {
    ////Part 1: test reservoir
    
    //Test_Grid<2>(); //a visual test, not included in general tests

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

    //uncomment those whey they don't rely on pre-saved data
    //Test_Mesh_Loader_Multiple<TriangleMesh<3>>(); // Can also verify by opening copy-mesh.obj
    //Test_Mesh_Loader_Single<TriangleMesh<3>>();

    /*Test_Center<float, 2>(Vector2i(328, 92));
    Test_Center<float, 3>(Vector3i(36, 28, 478));
    Test_Center<double, 2>(Vector2i(74, 22));
    Test_Center<double, 3>(Vector3i(381, 7, 32));*/

    ////Part 2: test sover

    Test_Sparse_Matrix();
    Test_CG_Memory_Safe();

    //Note: dense solver and sparse solver tests are in _tests_dec_system
}