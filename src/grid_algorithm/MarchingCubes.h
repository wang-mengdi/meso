//////////////////////////////////////////////////////////////////////////
// Marching Cubes Algorithm
// Copyright (c) (2022-), Bo Zhu, Yunquan Gu
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "Grid.h"
#include "Field.h"
#include "Mesh.h"

namespace Meso {

	__device__ static const int edge_table[256] = {
	0x0  , 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c,
	0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00,
	0x190, 0x99 , 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c,
	0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90,
	0x230, 0x339, 0x33 , 0x13a, 0x636, 0x73f, 0x435, 0x53c,
	0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30,
	0x3a0, 0x2a9, 0x1a3, 0xaa , 0x7a6, 0x6af, 0x5a5, 0x4ac,
	0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0,
	0x460, 0x569, 0x663, 0x76a, 0x66 , 0x16f, 0x265, 0x36c,
	0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60,
	0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0xff , 0x3f5, 0x2fc,
	0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0,
	0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x55 , 0x15c,
	0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950,
	0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0xcc ,
	0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0,
	0x8c0, 0x9c9, 0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc,
	0xcc , 0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9, 0x7c0,
	0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c,
	0x15c, 0x55 , 0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650,
	0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc,
	0x2fc, 0x3f5, 0xff , 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0,
	0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65, 0xc6c,
	0x36c, 0x265, 0x16f, 0x66 , 0x76a, 0x663, 0x569, 0x460,
	0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac,
	0x4ac, 0x5a5, 0x6af, 0x7a6, 0xaa , 0x1a3, 0x2a9, 0x3a0,
	0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c,
	0x53c, 0x435, 0x73f, 0x636, 0x13a, 0x33 , 0x339, 0x230,
	0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c,
	0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x99 , 0x190,
	0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c,
	0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x0
	};

	__device__ static const int triangle_table[256][16] =
	{ {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 8, 3, 9, 8, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 8, 3, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{9, 2, 10, 0, 2, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{2, 8, 3, 2, 10, 8, 10, 9, 8, -1, -1, -1, -1, -1, -1, -1},
	{3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 11, 2, 8, 11, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 9, 0, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 11, 2, 1, 9, 11, 9, 8, 11, -1, -1, -1, -1, -1, -1, -1},
	{3, 10, 1, 11, 10, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 10, 1, 0, 8, 10, 8, 11, 10, -1, -1, -1, -1, -1, -1, -1},
	{3, 9, 0, 3, 11, 9, 11, 10, 9, -1, -1, -1, -1, -1, -1, -1},
	{9, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{4, 3, 0, 7, 3, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 1, 9, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{4, 1, 9, 4, 7, 1, 7, 3, 1, -1, -1, -1, -1, -1, -1, -1},
	{1, 2, 10, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{3, 4, 7, 3, 0, 4, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1},
	{9, 2, 10, 9, 0, 2, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
	{2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, -1, -1, -1, -1},
	{8, 4, 7, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{11, 4, 7, 11, 2, 4, 2, 0, 4, -1, -1, -1, -1, -1, -1, -1},
	{9, 0, 1, 8, 4, 7, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
	{4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, -1, -1, -1, -1},
	{3, 10, 1, 3, 11, 10, 7, 8, 4, -1, -1, -1, -1, -1, -1, -1},
	{1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, -1, -1, -1, -1},
	{4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, -1, -1, -1, -1},
	{4, 7, 11, 4, 11, 9, 9, 11, 10, -1, -1, -1, -1, -1, -1, -1},
	{9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{9, 5, 4, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 5, 4, 1, 5, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{8, 5, 4, 8, 3, 5, 3, 1, 5, -1, -1, -1, -1, -1, -1, -1},
	{1, 2, 10, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{3, 0, 8, 1, 2, 10, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
	{5, 2, 10, 5, 4, 2, 4, 0, 2, -1, -1, -1, -1, -1, -1, -1},
	{2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, -1, -1, -1, -1},
	{9, 5, 4, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 11, 2, 0, 8, 11, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
	{0, 5, 4, 0, 1, 5, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
	{2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, -1, -1, -1, -1},
	{10, 3, 11, 10, 1, 3, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1},
	{4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10, -1, -1, -1, -1},
	{5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, -1, -1, -1, -1},
	{5, 4, 8, 5, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1},
	{9, 7, 8, 5, 7, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{9, 3, 0, 9, 5, 3, 5, 7, 3, -1, -1, -1, -1, -1, -1, -1},
	{0, 7, 8, 0, 1, 7, 1, 5, 7, -1, -1, -1, -1, -1, -1, -1},
	{1, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{9, 7, 8, 9, 5, 7, 10, 1, 2, -1, -1, -1, -1, -1, -1, -1},
	{10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, -1, -1, -1, -1},
	{8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, -1, -1, -1, -1},
	{2, 10, 5, 2, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1},
	{7, 9, 5, 7, 8, 9, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1},
	{9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, -1, -1, -1, -1},
	{2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7, -1, -1, -1, -1},
	{11, 2, 1, 11, 1, 7, 7, 1, 5, -1, -1, -1, -1, -1, -1, -1},
	{9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11, -1, -1, -1, -1},
	{5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0, -1},
	{11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0, -1},
	{11, 10, 5, 7, 11, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 8, 3, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{9, 0, 1, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 8, 3, 1, 9, 8, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
	{1, 6, 5, 2, 6, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 6, 5, 1, 2, 6, 3, 0, 8, -1, -1, -1, -1, -1, -1, -1},
	{9, 6, 5, 9, 0, 6, 0, 2, 6, -1, -1, -1, -1, -1, -1, -1},
	{5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, -1, -1, -1, -1},
	{2, 3, 11, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{11, 0, 8, 11, 2, 0, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
	{0, 1, 9, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
	{5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11, -1, -1, -1, -1},
	{6, 3, 11, 6, 5, 3, 5, 1, 3, -1, -1, -1, -1, -1, -1, -1},
	{0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, -1, -1, -1, -1},
	{3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, -1, -1, -1, -1},
	{6, 5, 9, 6, 9, 11, 11, 9, 8, -1, -1, -1, -1, -1, -1, -1},
	{5, 10, 6, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{4, 3, 0, 4, 7, 3, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1},
	{1, 9, 0, 5, 10, 6, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
	{10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, -1, -1, -1, -1},
	{6, 1, 2, 6, 5, 1, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1},
	{1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, -1, -1, -1, -1},
	{8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, -1, -1, -1, -1},
	{7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9, -1},
	{3, 11, 2, 7, 8, 4, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
	{5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11, -1, -1, -1, -1},
	{0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1},
	{9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6, -1},
	{8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6, -1, -1, -1, -1},
	{5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11, -1},
	{0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7, -1},
	{6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9, -1, -1, -1, -1},
	{10, 4, 9, 6, 4, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{4, 10, 6, 4, 9, 10, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1},
	{10, 0, 1, 10, 6, 0, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1},
	{8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10, -1, -1, -1, -1},
	{1, 4, 9, 1, 2, 4, 2, 6, 4, -1, -1, -1, -1, -1, -1, -1},
	{3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4, -1, -1, -1, -1},
	{0, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{8, 3, 2, 8, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1},
	{10, 4, 9, 10, 6, 4, 11, 2, 3, -1, -1, -1, -1, -1, -1, -1},
	{0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6, -1, -1, -1, -1},
	{3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10, -1, -1, -1, -1},
	{6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1, -1},
	{9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3, -1, -1, -1, -1},
	{8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1, -1},
	{3, 11, 6, 3, 6, 0, 0, 6, 4, -1, -1, -1, -1, -1, -1, -1},
	{6, 4, 8, 11, 6, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{7, 10, 6, 7, 8, 10, 8, 9, 10, -1, -1, -1, -1, -1, -1, -1},
	{0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10, -1, -1, -1, -1},
	{10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0, -1, -1, -1, -1},
	{10, 6, 7, 10, 7, 1, 1, 7, 3, -1, -1, -1, -1, -1, -1, -1},
	{1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, -1, -1, -1, -1},
	{2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9, -1},
	{7, 8, 0, 7, 0, 6, 6, 0, 2, -1, -1, -1, -1, -1, -1, -1},
	{7, 3, 2, 6, 7, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7, -1, -1, -1, -1},
	{2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7, -1},
	{1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11, -1},
	{11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1, -1, -1, -1, -1},
	{8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6, -1},
	{0, 9, 1, 11, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0, -1, -1, -1, -1},
	{7, 11, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{3, 0, 8, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 1, 9, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{8, 1, 9, 8, 3, 1, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
	{10, 1, 2, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 2, 10, 3, 0, 8, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
	{2, 9, 0, 2, 10, 9, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
	{6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8, -1, -1, -1, -1},
	{7, 2, 3, 6, 2, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{7, 0, 8, 7, 6, 0, 6, 2, 0, -1, -1, -1, -1, -1, -1, -1},
	{2, 7, 6, 2, 3, 7, 0, 1, 9, -1, -1, -1, -1, -1, -1, -1},
	{1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6, -1, -1, -1, -1},
	{10, 7, 6, 10, 1, 7, 1, 3, 7, -1, -1, -1, -1, -1, -1, -1},
	{10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8, -1, -1, -1, -1},
	{0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7, -1, -1, -1, -1},
	{7, 6, 10, 7, 10, 8, 8, 10, 9, -1, -1, -1, -1, -1, -1, -1},
	{6, 8, 4, 11, 8, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{3, 6, 11, 3, 0, 6, 0, 4, 6, -1, -1, -1, -1, -1, -1, -1},
	{8, 6, 11, 8, 4, 6, 9, 0, 1, -1, -1, -1, -1, -1, -1, -1},
	{9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6, -1, -1, -1, -1},
	{6, 8, 4, 6, 11, 8, 2, 10, 1, -1, -1, -1, -1, -1, -1, -1},
	{1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6, -1, -1, -1, -1},
	{4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9, -1, -1, -1, -1},
	{10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3, -1},
	{8, 2, 3, 8, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1},
	{0, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8, -1, -1, -1, -1},
	{1, 9, 4, 1, 4, 2, 2, 4, 6, -1, -1, -1, -1, -1, -1, -1},
	{8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1, -1, -1, -1, -1},
	{10, 1, 0, 10, 0, 6, 6, 0, 4, -1, -1, -1, -1, -1, -1, -1},
	{4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3, -1},
	{10, 9, 4, 6, 10, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{4, 9, 5, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 8, 3, 4, 9, 5, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
	{5, 0, 1, 5, 4, 0, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
	{11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5, -1, -1, -1, -1},
	{9, 5, 4, 10, 1, 2, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
	{6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5, -1, -1, -1, -1},
	{7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2, -1, -1, -1, -1},
	{3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6, -1},
	{7, 2, 3, 7, 6, 2, 5, 4, 9, -1, -1, -1, -1, -1, -1, -1},
	{9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7, -1, -1, -1, -1},
	{3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0, -1, -1, -1, -1},
	{6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8, -1},
	{9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7, -1, -1, -1, -1},
	{1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4, -1},
	{4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10, -1},
	{7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10, -1, -1, -1, -1},
	{6, 9, 5, 6, 11, 9, 11, 8, 9, -1, -1, -1, -1, -1, -1, -1},
	{3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5, -1, -1, -1, -1},
	{0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11, -1, -1, -1, -1},
	{6, 11, 3, 6, 3, 5, 5, 3, 1, -1, -1, -1, -1, -1, -1, -1},
	{1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6, -1, -1, -1, -1},
	{0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10, -1},
	{11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5, -1},
	{6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3, -1, -1, -1, -1},
	{5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2, -1, -1, -1, -1},
	{9, 5, 6, 9, 6, 0, 0, 6, 2, -1, -1, -1, -1, -1, -1, -1},
	{1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8, -1},
	{1, 5, 6, 2, 1, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6, -1},
	{10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0, -1, -1, -1, -1},
	{0, 3, 8, 5, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{10, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{11, 5, 10, 7, 5, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{11, 5, 10, 11, 7, 5, 8, 3, 0, -1, -1, -1, -1, -1, -1, -1},
	{5, 11, 7, 5, 10, 11, 1, 9, 0, -1, -1, -1, -1, -1, -1, -1},
	{10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1, -1, -1, -1, -1},
	{11, 1, 2, 11, 7, 1, 7, 5, 1, -1, -1, -1, -1, -1, -1, -1},
	{0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11, -1, -1, -1, -1},
	{9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7, -1, -1, -1, -1},
	{7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2, -1},
	{2, 5, 10, 2, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1},
	{8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5, -1, -1, -1, -1},
	{9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2, -1, -1, -1, -1},
	{9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2, -1},
	{1, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 8, 7, 0, 7, 1, 1, 7, 5, -1, -1, -1, -1, -1, -1, -1},
	{9, 0, 3, 9, 3, 5, 5, 3, 7, -1, -1, -1, -1, -1, -1, -1},
	{9, 8, 7, 5, 9, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{5, 8, 4, 5, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1},
	{5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0, -1, -1, -1, -1},
	{0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5, -1, -1, -1, -1},
	{10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4, -1},
	{2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8, -1, -1, -1, -1},
	{0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11, -1},
	{0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5, -1},
	{9, 4, 5, 2, 11, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4, -1, -1, -1, -1},
	{5, 10, 2, 5, 2, 4, 4, 2, 0, -1, -1, -1, -1, -1, -1, -1},
	{3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9, -1},
	{5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2, -1, -1, -1, -1},
	{8, 4, 5, 8, 5, 3, 3, 5, 1, -1, -1, -1, -1, -1, -1, -1},
	{0, 4, 5, 1, 0, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5, -1, -1, -1, -1},
	{9, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{4, 11, 7, 4, 9, 11, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1},
	{0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11, -1, -1, -1, -1},
	{1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11, -1, -1, -1, -1},
	{3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4, -1},
	{4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2, -1, -1, -1, -1},
	{9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3, -1},
	{11, 7, 4, 11, 4, 2, 2, 4, 0, -1, -1, -1, -1, -1, -1, -1},
	{11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4, -1, -1, -1, -1},
	{2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9, -1, -1, -1, -1},
	{9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7, -1},
	{3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10, -1},
	{1, 10, 2, 8, 7, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{4, 9, 1, 4, 1, 7, 7, 1, 3, -1, -1, -1, -1, -1, -1, -1},
	{4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1, -1, -1, -1, -1},
	{4, 0, 3, 7, 4, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{4, 8, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{9, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{3, 0, 9, 3, 9, 11, 11, 9, 10, -1, -1, -1, -1, -1, -1, -1},
	{0, 1, 10, 0, 10, 8, 8, 10, 11, -1, -1, -1, -1, -1, -1, -1},
	{3, 1, 10, 11, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 2, 11, 1, 11, 9, 9, 11, 8, -1, -1, -1, -1, -1, -1, -1},
	{3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9, -1, -1, -1, -1},
	{0, 2, 11, 8, 0, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{3, 2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{2, 3, 8, 2, 8, 10, 10, 8, 9, -1, -1, -1, -1, -1, -1, -1},
	{9, 10, 2, 0, 9, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8, -1, -1, -1, -1},
	{1, 10, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 3, 8, 9, 1, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 9, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1} };

	template<class T, int d>
	__host__ unsigned int Get_Cell_Type(const Field<T, d>& field, const Vector<int, d> cell_index, const real iso_value) {
		Typedef_VectorD(d);
		unsigned int cell_type = 0;
		if (field.Get(cell_index + VectorDi(0, 0, 0)) < iso_value) cell_type |= 1;
		if (field.Get(cell_index + VectorDi(1, 0, 0)) < iso_value) cell_type |= 2;
		if (field.Get(cell_index + VectorDi(1, 0, 1)) < iso_value) cell_type |= 4;
		if (field.Get(cell_index + VectorDi(0, 0, 1)) < iso_value) cell_type |= 8;
		if (field.Get(cell_index + VectorDi(0, 1, 0)) < iso_value) cell_type |= 16;
		if (field.Get(cell_index + VectorDi(1, 1, 0)) < iso_value) cell_type |= 32;
		if (field.Get(cell_index + VectorDi(1, 1, 1)) < iso_value) cell_type |= 64;
		if (field.Get(cell_index + VectorDi(0, 1, 1)) < iso_value) cell_type |= 128;
		return cell_type;
	}

	__host__ __device__ bool Has_Edge(const int edge_type, const int edge_index) { return (edge_type & (1 << edge_index)) == 0 ? false : true; }

	template<int d>
	__host__ __device__ void Get_Edge_Vertex_Index(const Vector<int, d>& cell_index, int edge_index, Vector<int, d>& v1, Vector<int, d>& v2)
	{
		Typedef_VectorD(d);
		switch (edge_index) {
		case 0:v1 = VectorDi(cell_index[0], cell_index[1], cell_index[2]);v2 = VectorDi(cell_index[0] + 1, cell_index[1], cell_index[2]);break;
		case 1:v1 = VectorDi(cell_index[0] + 1, cell_index[1], cell_index[2]);v2 = VectorDi(cell_index[0] + 1, cell_index[1], cell_index[2] + 1);break;
		case 2:v1 = VectorDi(cell_index[0] + 1, cell_index[1], cell_index[2]);v2 = VectorDi(cell_index[0] + 1, cell_index[1], cell_index[2] + 1);break;
		case 3:v1 = VectorDi(cell_index[0], cell_index[1], cell_index[2]);v2 = VectorDi(cell_index[0], cell_index[1], cell_index[2] + 1);break;
		case 4:v1 = VectorDi(cell_index[0], cell_index[1] + 1, cell_index[2]);v2 = VectorDi(cell_index[0] + 1, cell_index[1] + 1, cell_index[2]);break;
		case 5:v1 = VectorDi(cell_index[0] + 1, cell_index[1] + 1, cell_index[2]);v2 = VectorDi(cell_index[0] + 1, cell_index[1] + 1, cell_index[2] + 1);break;
		case 6:v1 = VectorDi(cell_index[0], cell_index[1] + 1, cell_index[2] + 1);v2 = VectorDi(cell_index[0] + 1, cell_index[1] + 1, cell_index[2] + 1);break;
		case 7:v1 = VectorDi(cell_index[0], cell_index[1] + 1, cell_index[2]);v2 = VectorDi(cell_index[0], cell_index[1] + 1, cell_index[2] + 1);break;
		case 8:v1 = VectorDi(cell_index[0], cell_index[1], cell_index[2]);v2 = VectorDi(cell_index[0], cell_index[1] + 1, cell_index[2]);break;
		case 9:v1 = VectorDi(cell_index[0] + 1, cell_index[1], cell_index[2]);v2 = VectorDi(cell_index[0] + 1, cell_index[1] + 1, cell_index[2]);break;
		case 10:v1 = VectorDi(cell_index[0] + 1, cell_index[1], cell_index[2] + 1);v2 = VectorDi(cell_index[0] + 1, cell_index[1] + 1, cell_index[2] + 1);break;
		case 11:v1 = VectorDi(cell_index[0], cell_index[1], cell_index[2] + 1);v2 = VectorDi(cell_index[0], cell_index[1] + 1, cell_index[2] + 1);break;
		default: break;
		}
	}

	template<int d>
	__host__ int& Index_On_Edge(const Vector<int, d>& cell_index, const int edge_index, Field<int, d>* v_idx_on_edge)
	{
		Typedef_VectorD(d);
		switch (edge_index) {
			////axis-x edges
		case 0:return v_idx_on_edge[0](VectorDi(cell_index[0], cell_index[1], cell_index[2]));
		case 2:return v_idx_on_edge[0](VectorDi(cell_index[0], cell_index[1], cell_index[2] + 1));
		case 4:return v_idx_on_edge[0](VectorDi(cell_index[0], cell_index[1] + 1, cell_index[2]));
		case 6:return v_idx_on_edge[0](VectorDi(cell_index[0], cell_index[1] + 1, cell_index[2] + 1));
			////axis-y edges
		case 8:return v_idx_on_edge[1](VectorDi(cell_index[0], cell_index[1], cell_index[2]));
		case 9:return v_idx_on_edge[1](VectorDi(cell_index[0] + 1, cell_index[1], cell_index[2]));
		case 11:return v_idx_on_edge[1](VectorDi(cell_index[0], cell_index[1], cell_index[2] + 1));
		case 10:return v_idx_on_edge[1](VectorDi(cell_index[0] + 1, cell_index[1], cell_index[2] + 1));
			////axis-z edges
		case 3:return v_idx_on_edge[2](VectorDi(cell_index[0], cell_index[1], cell_index[2]));
		case 1:return v_idx_on_edge[2](VectorDi(cell_index[0] + 1, cell_index[1], cell_index[2]));
		case 7:return v_idx_on_edge[2](VectorDi(cell_index[0], cell_index[1] + 1, cell_index[2]));
		case 5:return v_idx_on_edge[2](VectorDi(cell_index[0] + 1, cell_index[1] + 1, cell_index[2]));
		default: break;

		}
	}

	template<class T, int d>
	void Marching_Cubes_CPU(TriangleMesh<d> &_mesh, const Field<T, d>& field, const real _contour_val = 0.) {
		Typedef_VectorD(d); Typedef_VectorEi(d);
		
		const Grid<d>& grid = field.grid;
		const VectorDi cell_counts = field.grid.counts - VectorDi::Ones();

		Field<int, d> v_idx_on_edge[3]; for (int i = 0;i < 3;i++) v_idx_on_edge[i].Init(Grid<d>(grid.counts - VectorDi::Unit(i), grid.dx), -1);

		Array<VectorD>& vertices = *_mesh->vertices; vertices.clear();
		Array<VectorEi>& elements = _mesh->elements; elements.clear();

		grid.Exec_Nodes([&](const VectorDi cell_index) {
			if ((cell_index - cell_counts).maxCoeff() == 0) return;
			unsigned int cell_type = Get_Cell_Type<T, d>(field, cell_index, _contour_val);
			int edge_type = edge_table[cell_type];

			// go through all edges of one cell, then add vertex
			for (size_t ei = 0; ei < 12; ei++) {
				if (Has_Edge(edge_type, ei) && (Index_On_Edge<d>(cell_index, ei, v_idx_on_edge) == -1)) {
					VectorDi v1, v2; Get_Edge_Vertex_Index<d>(cell_index, ei, v1, v2);
					T alpha = (_contour_val - field.Get(v1)) / (field.Get(v2) - field.Get(v1));
					VectorD pos = (1 - alpha) * grid.Position(v1) + alpha * grid.Position(v2);
					vertices.push_back(pos); Index_On_Edge<d>(cell_index, ei, v_idx_on_edge) = (int)vertices.size() - 1;
				}

			}
			for (int ti = 0;triangle_table[cell_type][ti] != -1;ti += 3) {
				VectorDi t;
				t[0] = Index_On_Edge<d>(cell_index, triangle_table[cell_type][ti + 0], v_idx_on_edge);
				t[1] = Index_On_Edge<d>(cell_index, triangle_table[cell_type][ti + 1], v_idx_on_edge);
				t[2] = Index_On_Edge<d>(cell_index, triangle_table[cell_type][ti + 2], v_idx_on_edge);
				elements.push_back(t);
			}
			});

	}

	template<class T, int d>
	__device__ unsigned int Get_Cell_Type(const Grid<d> field_grid, const T* field_data, const Vector<int, d> cell_index, const real iso_value) {
		Typedef_VectorD(d);
		unsigned int cell_type = 0;
		if (field_data[field_grid.Index(cell_index + VectorDi(0, 0, 0))] < iso_value) cell_type |= 1;
		if (field_data[field_grid.Index(cell_index + VectorDi(1, 0, 0))] < iso_value) cell_type |= 2;
		if (field_data[field_grid.Index(cell_index + VectorDi(1, 0, 1))] < iso_value) cell_type |= 4;
		if (field_data[field_grid.Index(cell_index + VectorDi(0, 0, 1))] < iso_value) cell_type |= 8;
		if (field_data[field_grid.Index(cell_index + VectorDi(0, 1, 0))] < iso_value) cell_type |= 16;
		if (field_data[field_grid.Index(cell_index + VectorDi(1, 1, 0))] < iso_value) cell_type |= 32;
		if (field_data[field_grid.Index(cell_index + VectorDi(1, 1, 1))] < iso_value) cell_type |= 64;
		if (field_data[field_grid.Index(cell_index + VectorDi(0, 1, 1))] < iso_value) cell_type |= 128;
		return cell_type;
	}

	template<int d>
	__device__ int& Index_On_Edge(
		const Vector<int, d> cell_index,
		const int edge_index,
		const Grid<d>* edge_grid,
		int* v_idx_on_edge_x, int* v_idx_on_edge_y, int* v_idx_on_edge_z)
	{
		Typedef_VectorD(d);
		switch (edge_index) {
			////axis-x edges
		case 0: return v_idx_on_edge_x[edge_grid[0].Index(VectorDi(cell_index[0], cell_index[1], cell_index[2]))]; break;
		case 2: return v_idx_on_edge_x[edge_grid[0].Index(VectorDi(cell_index[0], cell_index[1], cell_index[2] + 1))]; break;
		case 4: return v_idx_on_edge_x[edge_grid[0].Index(VectorDi(cell_index[0], cell_index[1] + 1, cell_index[2]))]; break;
		case 6: return v_idx_on_edge_x[edge_grid[0].Index(VectorDi(cell_index[0], cell_index[1] + 1, cell_index[2] + 1))]; break;
			////axis-y edges
		case 8:  return v_idx_on_edge_y[edge_grid[1].Index(VectorDi(cell_index[0], cell_index[1], cell_index[2]))]; break;
		case 9:  return v_idx_on_edge_y[edge_grid[1].Index(VectorDi(cell_index[0] + 1, cell_index[1], cell_index[2]))]; break;
		case 11: return v_idx_on_edge_y[edge_grid[1].Index(VectorDi(cell_index[0], cell_index[1], cell_index[2] + 1))]; break;
		case 10: return v_idx_on_edge_y[edge_grid[1].Index(VectorDi(cell_index[0] + 1, cell_index[1], cell_index[2] + 1))]; break;
			////axis-z edges
		case 3: return v_idx_on_edge_z[edge_grid[2].Index(VectorDi(cell_index[0], cell_index[1], cell_index[2]))]; break;
		case 1: return v_idx_on_edge_z[edge_grid[2].Index(VectorDi(cell_index[0] + 1, cell_index[1], cell_index[2]))]; break;
		case 7: return v_idx_on_edge_z[edge_grid[2].Index(VectorDi(cell_index[0], cell_index[1] + 1, cell_index[2]))]; break;
		case 5: return v_idx_on_edge_z[edge_grid[2].Index(VectorDi(cell_index[0] + 1, cell_index[1] + 1, cell_index[2]))]; break;
		default: break;

		}
	}

	template<class T, int d>
	__global__ static void Mesh_Count(
		const Vector<int, d> cell_counts,
		const Grid<d> field_grid,
		const T* field_data,
		const real _contour_val,
		const Grid<d>* edge_grid,
		int* mesh_count, int* v_idx_on_edge_x, int* v_idx_on_edge_y, int* v_idx_on_edge_z
	) {
		Typedef_VectorD(d); Vector<int, d> cell_index = GPUFunc::Thread_Coord<d>(blockIdx, threadIdx);
		if ((cell_index - cell_counts).maxCoeff() == 0) return;
		unsigned int cell_type = Get_Cell_Type<T, d>(field_grid, field_data, cell_index, _contour_val);
		int edge_type = edge_table[cell_type];

		// Reserve space for mesh
		for (int ti = 0;triangle_table[cell_type][ti * 3] != -1;ti += 1) mesh_count[field_grid.Index(cell_index)] += 1;
		// Go through all edges of the cell, and label existing point
		for (int ei = 0; ei < 12; ei++) {
			if (Has_Edge(edge_type, ei))
				Index_On_Edge<d>(cell_index, ei, edge_grid, v_idx_on_edge_x, v_idx_on_edge_y, v_idx_on_edge_z) = 1;
		}

	}

	template<class T, int d>
	__global__ static void Gen_Vertice(
		const Vector<int, d> cell_counts,
		const Grid<d> field_grid,
		const T* field_data,
		const real _contour_val,
		const Grid<d> edge_grid,
		int axis,
		int* v_idx_on_edge, Vector<real, d>* vertices)
	{
		Typedef_VectorD(d); Vector<int, d> cell_index = GPUFunc::Thread_Coord<d>(blockIdx, threadIdx);
		int index = v_idx_on_edge[edge_grid.Index(cell_index)];
		if (index >= 0) {
			VectorDi v1 = cell_index; VectorDi v2 = v1 + VectorDi::Unit(axis);
			T alpha = (_contour_val - field_data[field_grid.Index(v1)]) / (field_data[field_grid.Index(v2)] - field_data[field_grid.Index(v1)]);
			VectorD pos = (1 - alpha) * field_grid.Position(v1) + alpha * field_grid.Position(v2);
			vertices[index] = pos;
		}
	}
	
	template<class T, int d>
	__global__ void Generate_Mesh_Kernel(
		const Vector<int, d> cell_counts,
		const Grid<d> field_grid,
		const T* field_data,
		const real _contour_val,
		const int* window,
		const Grid<d>* edge_grid,
		int* v_idx_on_edge_x, int* v_idx_on_edge_y, int* v_idx_on_edge_z,
		Vector<int, d>* meshes
	) {
		Typedef_VectorD(d); Vector<int, d> cell_index = GPUFunc::Thread_Coord<d>(blockIdx, threadIdx);
		if ((cell_index - cell_counts).maxCoeff() == 0) return;
		unsigned int cell_type = Get_Cell_Type<T, d>(field_grid, field_data, cell_index, _contour_val);
		int edge_type = edge_table[cell_type];

		int start = window[field_grid.Index(cell_index)];
		for (int ti = 0;triangle_table[cell_type][ti * 3] != -1;ti += 1) {
			meshes[start + ti][0] = Index_On_Edge<d>(cell_index, triangle_table[cell_type][ti * 3 + 0], edge_grid, v_idx_on_edge_x, v_idx_on_edge_y, v_idx_on_edge_z);
			meshes[start + ti][1] = Index_On_Edge<d>(cell_index, triangle_table[cell_type][ti * 3 + 1], edge_grid, v_idx_on_edge_x, v_idx_on_edge_y, v_idx_on_edge_z);
			meshes[start + ti][2] = Index_On_Edge<d>(cell_index, triangle_table[cell_type][ti * 3 + 2], edge_grid, v_idx_on_edge_x, v_idx_on_edge_y, v_idx_on_edge_z);
		}
	}

	template<class T, int d>
	void Marching_Cubes_GPU(
		TriangleMesh<d>& _mesh,
		const Field<T, d, DEVICE>& field,
		const real _contour_val = 0.
	) {
		Typedef_VectorD(d); Typedef_VectorEi(d);

		// 0. Init
		const Grid<d>& grid = field.grid; _mesh = (_mesh == nullptr) ? std::make_shared<TriangleMesh<3> >() : _mesh;
		const VectorDi cell_counts = grid.counts - VectorDi::Ones();

		ArrayDv<int> mesh_count(grid.DoF() + 1, 0);

		Array<Grid<d>> edge_grid(3); for (int i = 0;i < 3;i++) edge_grid[i] = Grid<d>(grid.counts - VectorDi::Unit(i), grid.dx);
		ArrayDv<Grid<d>> edge_grid_dv = edge_grid;

		FieldDv<int, d> v_idx_on_edge[3]; Field<int, d> v_idx_on_edge_host[3];
		for (int i = 0;i < 3;i++) v_idx_on_edge[i].Init(edge_grid[i], -1);

		// 1. Reserve memory for mesh
		grid.Exec_Kernel(
			&Mesh_Count<T, d>,
			cell_counts, grid, field.Data_Ptr(), _contour_val, ArrayFunc::Data<Grid<d>, DEVICE>(edge_grid_dv), // const data
			ArrayFunc::Data<int, DEVICE>(mesh_count), v_idx_on_edge[0].Data_Ptr(), v_idx_on_edge[1].Data_Ptr(), v_idx_on_edge[2].Data_Ptr()
		);

		// 2. Generate vertices
		size_t total = 0; 
#pragma omp parallel for
		for (size_t i = 0; i < d; i++) {
			size_t count = thrust::count(v_idx_on_edge[i].data->begin(), v_idx_on_edge[i].data->end(), 1);
			thrust::device_vector<int> seq(v_idx_on_edge[0].data->size()); thrust::sequence(thrust::device, seq.begin(), seq.end(), 0);
			thrust::device_vector<int> p_iter(count);
			thrust::copy_if(seq.begin(), seq.end(), v_idx_on_edge[i].data->begin(), p_iter.begin(), thrust::placeholders::_1 > 0);
			thrust::copy_n(thrust::make_counting_iterator(total), count, thrust::make_permutation_iterator(v_idx_on_edge[i].data->begin(), p_iter.begin()));
			total += count;
		}
		
		Array<VectorD, DEVICE> vertices(total);
#pragma omp parallel for
		for (size_t i = 0; i < d; i++) edge_grid[i].Exec_Kernel(
			&Gen_Vertice<T, d>,
			cell_counts, grid, field.Data_Ptr(), _contour_val, edge_grid_dv[i], i, // const data
			v_idx_on_edge[i].Data_Ptr(), ArrayFunc::Data<VectorD, DEVICE>(vertices)
		);

		// 3. Generate Meshes
		thrust::device_ptr<int> mesh_count_ptr = thrust::device_pointer_cast(mesh_count.data());
		thrust::exclusive_scan(mesh_count_ptr, mesh_count_ptr + mesh_count.size(), mesh_count_ptr);
		Array<VectorEi, DEVICE> elements(mesh_count[mesh_count.size() - 1]);

		grid.Exec_Kernel(
			&Generate_Mesh_Kernel<T, d>,
			cell_counts,
			grid, field.Data_Ptr(), _contour_val,
			ArrayFunc::Data<int, DEVICE>(mesh_count),
			ArrayFunc::Data<Grid<d>, DEVICE>(edge_grid_dv),
			v_idx_on_edge[0].Data_Ptr(), v_idx_on_edge[1].Data_Ptr(), v_idx_on_edge[2].Data_Ptr(),
			ArrayFunc::Data<VectorEi, DEVICE>(elements)
		);

		*_mesh->vertices = vertices;
		_mesh->elements = elements;
	}

	template<class T, int d, DataHolder side = HOST>
	void Marching_Cubes(
		TriangleMesh<d>& _mesh,
		const Field<T, d, side>& field,
		const real _contour_val = 0.
	) {
		if constexpr (side == HOST) {
			Marching_Cubes_CPU(_mesh, field, _contour_val);
		}
		else {
			Marching_Cubes_GPU(_mesh, field, _contour_val);
		}
	}
}