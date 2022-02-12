#include <iostream>
#include <cstdio>

#include "Common.h"
#include "Grid.h"
#include "Field.h"
#include <fmt/ranges.h>

template<int d>
void Test_Grid(void) {
    Typedef_VectorD(d);
    Grid<2> grid(Vector2i(9, 7));
    Field<int, 2> F(grid);
    grid.Exec_Each(
        [&](const VectorDi& cell) {

        }
    );
}

int main(){
    Assert(false, "alert");
    //Test_Grid<2>();
    return 0;
}