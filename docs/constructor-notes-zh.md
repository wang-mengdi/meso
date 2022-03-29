# 构造函数的小括号歧义
如下代码无法编译：
```c++
#include <iostream>
#include <typeinfo>
using namespace std;

class Grid{
    public:
    int size;
    Grid(const int _size=0):size(_size){}
};

template<int d>
class Field{
    public:
    Grid grid;
    Field(){
    }
    Field(const Grid &_grid){
        grid=_grid;
    }
};

void test(int sz){
    Field<2> af(Grid(sz));
    std::cout<<af.grid.size<<"\n";
}

int main(){
    test(3);
    return 0;
}
```

编译器提示这一行报错：
```c++
std::cout<<af.grid.size<<"\n";
```
提示`af`的类型是`Field<2>(Grid)`。其含义实际上是，编译器把`af`的定义：
```c++
Field<2> af(Grid(sz));
```
视作了一个函数声明，即一个名为`af`的函数，输入类型为`Grid(int)`，输出类型为`Field<2>`。有以下几种改动可以解决这一问题：
1. 用大括号调用构造函数：
```c++
Field<2> af{Grid<sz>};
```
2. 如果构造函数传入一个右值，编译器不会将其视作函数声明。例如：
```c++
Field<2> af(Grid<3>);
Field<2> af(std::move(Grid<sz>));
```
3. 如果`Field`类无模板参数`d`，也不会出现这一问题，可能是编译器的自动识别。
