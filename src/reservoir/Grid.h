//////////////////////////////////////////////////////////////////////////
// Basic grid data structure
// Copyright (c) (2018-), Bo Zhu, Zangyueyang Xian, Mengdi Wang
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once

template<int d>
class Grid {
public:
	static constexpr int block_size = (d == 2 ? 8 : 4);
};