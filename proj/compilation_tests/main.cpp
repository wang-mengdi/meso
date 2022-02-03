#include <iostream>
#include <cstdio>
#include "Common.h"
int main(){
    foo();
    int sum = 0;
    for (int i = 0; i < 100000000; i++) {
        sum += i * i;
    }
    printf("%d\n", sum);
    return 0;
}