# 第三方库
## fmt

目前nvcc编译器不支持最新版的fmt(9.1.0)。如将<fmt/ostream.h> include到.cu文件中 或者 在.cu文件中调用<fmt/color.h>中的fmt::print()会导致编译错误。在nvcc编译器更新以前，先将fmt的版本设置为8.1.1。
