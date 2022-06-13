//////////////////////////////////////////////////////////////////////////
// Mesh
// Copyright (c) (2022-) Bo zhu, Yunquan Gu, Mengdi Wang, Fan Feng
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include <Common.h>

namespace Meso {
	template<class T, int d> using VertexMatrix = Eigen::Matrix<T, Eigen::Dynamic, d, Eigen::RowMajor>; //each row is the coordinate of a triangle
	template<int e_d> using ElementMatrix = Eigen::Matrix<int, Eigen::Dynamic, e_d, Eigen::RowMajor>; //each row contains the indices of one element
}