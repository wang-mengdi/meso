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
            F(cell) = grid.Index(cell);
        }
    );
    for (int i = 0; i < grid.counts[0]; i++) {
        for (int j = 0; j < grid.counts[1]; j++) {
            fmt::print("{} ", F(Vector2i(i, j)));
        }
        fmt::print("\n");
    }
}

int main(){
    Test_Grid<2>();
    return 0;
}