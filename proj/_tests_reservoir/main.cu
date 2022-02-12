#include <iostream>
#include <cstdio>

#include "Common.h"
#include "Grid.h"
#include "Field.h"
#include <fmt/ranges.h>

void Test_Grid(void) {
    Grid<2> grid(Vector2i(9, 7));
    Field<int, 2> F(grid);

}

int main(){
    //fmt::print(fg(fmt::color::red), "wtfwtftwt fwtf");
    //std::cout << "\n";
    Info("some info: {}", Vector2(9.66545335,8.7456456456));
    Grid<2> a(Vector2i(2, 5));
    Info("a size: {}", a.block_size);
    Info("size of vector3i: {}", sizeof(Vector3i));

    Vector3i b(1, 2, 4), c(1, 2, 3);
    Warn("b==c: {}", b == c);
    return 0;
}