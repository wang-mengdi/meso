# 编程注意事项

- 使用`Json::Value()`解析字符串时，必须显式指定模板std::string，即`Json::Value<std::string>(...);`，或把默认值参数显式转化为std::string，即`Json::Value(j,"parameter",(std::string)"some_string");`，否则会报一个模板解析错误，大概是什么什么东西无法和const char*转换。
- `Interpolation`中用到了`FaceField`，而`FaceField`中的`Output_Vtk`使用了`Interpolation`，此处循环定义通过在`Interpolation.h`中对`FaceField`类的forward declare解决。若不解决，会报一个"FaceField不是模板"的错误。
- 目前无法使用Eigen的incomplete cholesky preconditioner，会报一个“无法访问private typedef”的错误。
- 目前我们无法在Eigen的3.4.0版本中使用Quaternion或者Matrix3初始化一个AngleAxis类。
- 目前OBJFunc使用tinyobjload的1.0.7版本（xmake repo版本原因），等xmake更新后再换成2.0.0的API
- 如果某个类的某个构造函数的参数中没有template，那么就不能在这个构造函数的定义中添加模板。此事出现过一次：DampedJacobiSmoother类曾经有一个带模板int d的构造函数，用于用PoissonMapping初始化，后来此构造函数被移除，但模板仍然留存，这就导致编译器报了一些无法识别的编译错误，错误提示无法匹配构造函数，但不知道应该如何匹配。
- 当调用parent class的template函数foo<T>的时候需要用parent.template foo<T>();
- 当在src创建新的library的时候，确保有一个.cu或.cpp file。否则没有.lib文件会被创建。
- 在`SparseSolverTests.h`中，为了支持:

    ```C++
    CholeskySparseSolver<T> solver(SparseMatrix_From_PoissonLike(grid, poisson_mapping));
    ```
    把`CholeskySparseSolver`和`SparseMatrixMapping`的构造函数都设置成了传入右值。右值是一个比较tricky的领域，暂时不列入标准，只是进行一些探索。
- 在函数参数里面放一个像`teplate<T> Array<T>`这样的东西，一般无法自动推导，因为`thrust::device_vector`或者`host_vector`都有两个模板参数，第一个是`T`，第二个是`allocator<T>`.