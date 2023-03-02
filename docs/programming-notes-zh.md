# 编程注意事项

- 使用`Json::Value()`解析字符串时，必须显式指定模板std::string，即`Json::Value<std::string>(...);`，或把默认值参数显式转化为std::string，即`Json::Value(j,"parameter",(std::string)"some_string");`，否则会报一个模板解析错误，大概是什么什么东西无法和const char*转换。
- `Interpolation`中用到了`FaceField`，而`FaceField`中的`Output_Vtk`使用了`Interpolation`，此处循环定义通过在`Interpolation.h`中对`FaceField`类的forward declare解决。若不解决，会报一个"FaceField不是模板"的错误。
- 目前无法使用Eigen的incomplete cholesky preconditioner，会报一个“无法访问private typedef”的错误。这个在simplex里面是通过自己魔改eigen实现的，因此为兼容考虑，在96515390b7d1bc36d527a76224b79c92a629a03e不再以public模式引用eigen包，而是在每一个xmake.lua里面非public地add_package一遍。
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
- 在函数参数里面放一个像`template<T> Array<T>`这样的东西，一般无法自动推导，因为`thrust::device_vector`或者`host_vector`都有两个模板参数，第一个是`T`，第二个是`allocator<T>`.
- 在求解泊松系统时，如果边界条件封闭（即周围一圈都是墙），那么意味着对应的线性系统有非平凡零空间（即，可以把流体压强加上任意常数$C$而不改变结果）。特别地，如果没有任何受力，那么`MultiGrid`最粗一层的直接求解会十分倾向于解出来一个错误的压强，甚至直接出现`NaN`.根据目前我们的观察，这种问题仅仅会在这一个特定的情况下出现，一旦流体“动起来”基本就没有问题了，如果调试遇到不必特别在意。
- 在global函数中，特定情况下不能使用VectorD，应显示写出`Vector<T,d>`，例如以下情况，不然数值会不变。具体原因不详。
  ```template<class T, int d>
	__global__ static void Covector_Stretching_Cell(const Grid<d> grid, T* result_data, const Vector<T, d>* inverse_flow_map_grad, 
		const Grid<d> gv0, const T* v0, const Grid<d> gv1, const T* v1, const Grid<d> gv2, const T* v2) {
		Typedef_VectorD(d);
		Vector<int, d> cell = GPUFunc::Thread_Coord<d>(blockIdx, threadIdx);
		VectorD pos = grid.Position(cell);
		Vector<T, d> grad = inverse_flow_map_grad[grid.Index(cell)];
		Vector<T, d> vel = IntpLinearPadding0::Face_Vector(gv0, v0, gv1, v1, gv2, v2, pos);
		result_data[grid.Index(cell)] = grad.dot(vel);
	}
    ```
- Add_Scalar 函数可以用于 Array\<Vector> 的计算，但是HOST和DEVICE会需要一些不同的处理，如直接用在使用HOST时n_a_i，则此函数并无效果。
``` VectorD a_i = a[i];
	VectorD n_a_i = -a_i;
	if constexpr (side == HOST) {
		ArrayFunc::Add_Scalar(temp_b, -a_i);
	}
	else {
		ArrayFunc::Add_Scalar(temp_b, n_a_i);
	}
```
