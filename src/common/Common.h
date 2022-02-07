//////////////////////////////////////////////////////////////////////////
// Common header
// Copyright (c) (2018-), Bo Zhu, Mengdi Wang
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once

#include "Eigen/Dense"
#include "Eigen/Geometry"
#include <fmt/core.h>

#include <iostream>

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


#ifdef USE_FLOAT
Declare_Eigen_Types(float, f)
#else
Declare_Eigen_Types(double, d)
#endif

#define Typedef_VectorD(d) \
using VectorD=Vector<real,d>; \
using VectorDi=Vector<int,d>

//// fmt part

template <typename... Args>
void Assert(const bool flg, const char* fmt = "", const Args &...args) {
    if (!flg) {
        std::cerr << "[Error]" << fmt::format(fmt, args...);
        exit(-1);
    }
}

template<typename ...Args>
void Info(const char* fmt, const Args&...args) {
    std::cout << "#     ";
    fmt::print(fmt, args...);
    std::cout << "\n";
}

void Info(const std::string& str);