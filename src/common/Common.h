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

#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "Eigen/Geometry"
#include <fmt/core.h>
#include <fmt/color.h>
#include <fmt/ostream.h>

#include <iostream>
#include <vector>

namespace Meso {
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

    Declare_Eigen_Vector_Types(int, i);
    Declare_Eigen_Vector_Types(float, f);
    Declare_Eigen_Vector_Types(double, d);

    using uchar = unsigned char;
    using ushort = unsigned short;
    template<class T, int d> using Vector = Eigen::Matrix<T, d, 1>;
    template<class T, int d> using Matrix = Eigen::Matrix<T, d, d>;

#define Typedef_VectorD(d) \
using VectorD=Vector<real,d>; \
using VectorDi=Vector<int,d>

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

    //// fmt part

    template<typename ...Args>
    void Info(const char* fmt, const Args&...args) {
        fmt::print("#     ");
        fmt::print(fmt, args...);
        fmt::print("\n");
    }
    void Info(const std::string& str);

    template<typename ...Args>
    void Warn(const char* fmt, const Args&...args) {
        fmt::print(fg(fmt::color::yellow), "#     ");
        fmt::print(fg(fmt::color::yellow), fmt, args...);
        fmt::print("\n");
    }
    void Warn(const std::string& str);

    template<typename ...Args>
    void Error(const char* fmt, const Args&...args) {
        fmt::print(fg(fmt::color::red), "#     ");
        fmt::print(fg(fmt::color::red), fmt, args...);
        fmt::print("\n");
    }
    void Error(const std::string& str);

    template<typename ...Args>
    void Pass(const char* fmt, const Args&...args) {
        fmt::print(fg(fmt::color::green), "#     ");
        fmt::print(fg(fmt::color::green), fmt, args...);
        fmt::print("\n");
    }
    void Pass(const std::string& str);

    template <typename... Args>
    void Assert(const bool flg, const char* fmt = "", const Args &...args) {
        if (!flg) {
            Error(fmt, args...);
            exit(-1);
        }
    }

    //////fmt adaptor for eigen vector
    //template <class T, int d> struct fmt::formatter<Vector<T, d> > {
    //    constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin()) {
    //        //https://fmt.dev/latest/api.html#udt
    //        auto it = ctx.begin(), end = ctx.end();
    //        if (it != end && *it != '}') throw format_error("invalid format");

    //        // Return an iterator past the end of the parsed range:
    //        return it;
    //    }

    //    // Formats the point p using the parsed format specification (presentation)
    //    // stored in this formatter.
    //    template <typename FormatContext>
    //    auto format(const Eigen::Matrix<T, d, 1>& vec, FormatContext& ctx) -> decltype(ctx.out()) {
    //        std::stringstream ss;
    //        ss << vec.transpose();
    //        // ctx.out() is an output iterator to write to.
    //        return format_to(
    //            ctx.out(),
    //            "{}",
    //            ss.str());
    //    }
    //};
}

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