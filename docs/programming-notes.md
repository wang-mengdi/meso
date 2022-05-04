# 编程注意事项

- 使用`Json::Value()`解析字符串时，必须显式指定模板std::string，即`Json::Value<std::string>(...);`，或把默认值参数显式转化为std::string，即`Json::Value(j,"parameter",(std::string)"some_string");`，否则会报一个模板解析错误，大概是什么什么东西无法和const char*转换。
- `Interpolation`中用到了`FaceField`，而`FaceField`中的`Output_Vtk`使用了`Interpolation`，此处循环定义通过在`Interpolation.h`中对`FaceField`类的forward declare解决。若不解决，会报一个"FaceField不是模板"的错误。
- 目前无法使用Eigen的incomplete cholesky preconditioner，会报一个“无法访问private typedef”的错误。
- 目前我们无法在Eigen的3.4.0版本中使用Quaternion或者Matrix3初始化一个AngleAxis类。