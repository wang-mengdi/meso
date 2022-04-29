//////////////////////////////////////////////////////////////////////////
// Common header
// Copyright (c) (2022-), Bo Zhu, Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once

#include <cuda_runtime.h>
#include "cublas_v2.h"
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <vector_functions.h>
#include "driver_types.h"

#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "Eigen/Geometry"
#include <fmt/core.h>
#include <fmt/color.h>
#include <fmt/ostream.h>
//fmt/ranges.h will override the format of Vector<T,d>
#include <fmt/ranges.h>
//ideally we don't want to use standard list, queue and array
#include <list>
#include <queue>
#include <array>
#include <iostream>
#include <vector>
#include <cmath>

#include <boost/filesystem.hpp>

namespace Meso {
    namespace bf = boost::filesystem;

    ////template operations
    template<bool condition, class T1, class T2> struct If { typedef T1 Type; };
    template<class T1, class T2> struct If<false, T1, T2> { typedef T2 Type; };

    ////Eigen part

#define Declare_Eigen_Types(my_type,t)  \
using real=my_type;                     \
using Vector1=Eigen::Matrix<real,1,1>;  \
using Vector2=Eigen::Vector2##t;        \
using Vector3=Eigen::Vector3##t;        \
using Vector4=Eigen::Vector4##t;        \
using VectorX=Eigen::VectorX##t;        \
using Matrix2=Eigen::Matrix2##t;        \
using Matrix3=Eigen::Matrix3##t;        \
using Matrix4=Eigen::Matrix4##t;        \
using MatrixX=Eigen::MatrixX##t;        \
using C=std::complex<real>;             \
using Quaternion=Eigen::Quaternion##t;  \
using AngleAxis=Eigen::AngleAxis##t;	\

    template<class T, int d> using Matrix = Eigen::Matrix<T, d, d>;
    template<class T> using SparseMatrix = Eigen::SparseMatrix<T, Eigen::RowMajor, int>;

#ifdef USE_FLOAT
    Declare_Eigen_Types(float, f)
#else
    Declare_Eigen_Types(double, d)
#endif

#define Declare_Eigen_Vector_Types(type,t)		\
using Vector1##t=Eigen::Matrix<type,1,1>;       \
using Vector2##t=Eigen::Vector2##t;             \
using Vector3##t=Eigen::Vector3##t;             \
using Vector4##t=Eigen::Vector4##t;             \
using VectorX##t=Eigen::VectorX##t;				

#define Declare_Eigen_Matrix_Types(type,t)		\
using Matrix1##t=Eigen::Matrix<type,1,1>;       \
using Matrix2##t=Eigen::Matrix2##t;             \
using Matrix3##t=Eigen::Matrix3##t;             \
using Matrix4##t=Eigen::Matrix4##t;             \
using MatrixX##t=Eigen::MatrixX##t;             

Declare_Eigen_Vector_Types(int, i);
Declare_Eigen_Vector_Types(float, f);
Declare_Eigen_Vector_Types(double, d);
Declare_Eigen_Matrix_Types(int, i)
Declare_Eigen_Matrix_Types(float, f)
Declare_Eigen_Matrix_Types(double, d)

using uchar = unsigned char;
using ushort = unsigned short;
template<class T, int d> using Vector = Eigen::Matrix<T, d, 1>;
template<class T, int d> using Matrix = Eigen::Matrix<T, d, d>;

#define Typedef_VectorD(d) \
using VectorD=Vector<real,d>; \
using VectorDi=Vector<int,d>; 
#define Typedef_MatrixD(d) \
using MatrixD=Matrix<real,d>;
#define Typedef_VectorEi(d) \
using VectorEi=Vector<int,d>;

// CUDA programming
enum DataHolder { HOST = 0, DEVICE };

////Container alias

//Array
template<class T, DataHolder side = DataHolder::HOST>
using Array = typename std::conditional<side == DataHolder::HOST, thrust::host_vector<T>, thrust::device_vector<T>>::type;
//template<class T> using Array = thrust::host_vector<T>;
template<class T> using ArrayDv = thrust::device_vector<T>;//device array
template<class T, DataHolder side = DataHolder::HOST> using ArrayPtr = std::shared_ptr<Array<T, side> >;
template<class T> using ArrayDvPtr = std::shared_ptr<ArrayDv<T> >;//device array ptr
////Array with fixed size
template<class T, int n> using ArrayF = std::array<T, n>;
constexpr int Pow(int x, int p) { return p == 1 ? x : x * Pow(x, p - 1); }
template<class T, int d> using ArrayF2P = ArrayF<T, Pow(2, d) >;
template<class T, int d> using ArrayF3P = ArrayF<T, Pow(3, d) >;
////Other containers
template<class T> using List = std::list<T>;
template<class T, class CMP = std::less<T> > using Heap = std::priority_queue<T, std::vector<T>, CMP>;

static const char* _cudaGetErrorEnum(cudaError_t error) {
    return cudaGetErrorName(error);
}
#define checkCudaErrors(val) check((val), #val, __FILE__, __LINE__)
    template <typename T>
    void check(T result, char const* const func, const char* const file,
        int const line) {
        if (result) {
            fprintf(stderr, "CUDA error at %s:%d code=%d(%s) \"%s\" \n", file, line,
                static_cast<unsigned int>(result), _cudaGetErrorEnum(result), func);
            exit(EXIT_FAILURE);
        }
    }

    void Check_Cuda_Memory(const std::string str = "");

    //// fmt part

    template<class... Args>
    void Info(const char* fmt_str, const Args&...args) {
        fmt::print("#     ");
        //auto fst = fmt::format_string<Args...>(fmt_str);
        //fmt::print(fst, (args)...);
        fmt::print(fmt_str, args...);
        fmt::print("\n");
    }
    void Info(const std::string& str);

    template<typename ...Args>
    void Warn(const char* fmt_str, const Args&...args) {
        fmt::print(fg(fmt::color::yellow), "#     ");
        fmt::print(fg(fmt::color::yellow), fmt_str, args...);
        fmt::print("\n");
    }
    void Warn(const std::string& str);

    template<typename ...Args>
    void Error(const char* fmt_str, const Args&...args) {
        fmt::print(fg(fmt::color::red), "#     ");
        fmt::print(fg(fmt::color::red), fmt_str, args...);
        fmt::print("\n");
    }
    void Error(const std::string& str);

    template<typename ...Args>
    void Pass(const char* fmt_str, const Args&...args) {
        fmt::print(fg(fmt::color::green), "#     ");
        fmt::print(fg(fmt::color::green), fmt_str, args...);
        fmt::print("\n");
    }
    void Pass(const std::string& str);

    template <typename... Args>
    void Assert(const bool flg, const char* fmt_str = "", const Args &...args) {
        if (!flg) {
            Error(fmt_str, args...);
            exit(-1);
        }
    }
}

//////fmt adaptor for eigen vector
//////not compatible with fmt/range.h
//template <class T, int d> 
//struct fmt::formatter<Meso::Vector<T, d> > {
//    constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin()) {
//        //https://fmt.dev/latest/api.html#udt
//        auto it = ctx.begin(), end = ctx.end();
//        if (it != end && *it != '}') throw format_error("invalid format");
//
//        // Return an iterator past the end of the parsed range:
//        return it;
//    }
//
//    // Formats the point p using the parsed format specification (presentation)
//    // stored in this formatter.
//    template <typename FormatContext>
//   // auto format(const Eigen::Matrix<T, d, 1>& vec, FormatContext& ctx) -> decltype(ctx.out()) {
//    auto format(const Meso::Vector<T, d>& vec, FormatContext& ctx) -> decltype(ctx.out()) {
//        std::stringstream ss;
//        ss << vec.transpose();
//        // ctx.out() is an output iterator to write to.
//        return format_to(
//            ctx.out(),
//            "{}",
//            ss.str());
//    }
//};

////fmt adaptor for eigen vector
//template <class T> struct fmt::formatter<Eigen::Matrix<T, 3, 3> > {
//    constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin()) {
//        //https://fmt.dev/latest/api.html#udt
//        auto it = ctx.begin(), end = ctx.end();
//        if (it != end && *it != '}') throw format_error("invalid format");
//
//        // Return an iterator past the end of the parsed range:
//        return it;
//    }
//
//    // Formats the point p using the parsed format specification (presentation)
//    // stored in this formatter.
//    template <typename FormatContext>
//    auto format(const Eigen::Matrix<T, 3, 3>& mat, FormatContext& ctx) -> decltype(ctx.out()) {
//        std::stringstream ss;
//        ss << mat;
//        // ctx.out() is an output iterator to write to.
//        return format_to(
//            ctx.out(),
//            "{}",
//            ss.str());
//    }
//};
//
//template<class T> struct fmt::formatter<Eigen::Quaternion<T>> {
//    constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin()) {
//        auto it = ctx.begin(), end = ctx.end();
//        if (it != end && *it != '}') throw format_error("invalid format");
//
//        // Return an iterator past the end of the parsed range:
//        return it;
//    }
//    template <typename FormatContext>
//    auto format(const Eigen::Quaternion<T>& quaternion, FormatContext& ctx) -> decltype(ctx.out()) {
//        std::stringstream ss;
//        ss << quaternion.coeffs().transpose();
//        // ctx.out() is an output iterator to write to.
//        return format_to(
//            ctx.out(),
//            "{}",
//            ss.str());
//    }
//};