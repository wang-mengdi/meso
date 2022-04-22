# 编程注意事项

- 使用`Json::Value()`解析字符串时，必须显式指定模板std::string，即`Json::Value<std::string>(...);`，或把默认值参数显式转化为std::string，即`Json::Value(j,"parameter",(std::string)"some_string");`，否则会报一个模板解析错误，大概是什么什么东西无法和const char*转换。
