#include <iostream>
#include <cstdio>
#include "Common.h"
#include "Grid.h"
#include <fmt/ranges.h>

int main(){
    Grid<2> a;
    Info("a size: {}", Grid<2>::block_size);
    return 0;
}