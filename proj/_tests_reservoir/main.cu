#include <iostream>
#include <cstdio>

#include "Common.h"
#include "Grid.h"
#include <fmt/ranges.h>

int main(){
    fmt::print(fg(fmt::color::red), "wtfwtftwt fwtf");
    std::cout << "\n";
    Grid<real, 2> a(Vector2i(2, 5));
    Info("a size: {}", a.block_size);
    Info("size of vector3i: {}", sizeof(Vector3i));

    Vector3i b(1, 2, 4), c(1, 2, 3);
    Warn("b==c: {}", b == c);
    return 0;
}