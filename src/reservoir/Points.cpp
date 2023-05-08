//////////////////////////////////////////////////////////////////////////
// Basic grid data representation
// Copyright (c) (2022-), Bo Zhu, Mengdi Wang, Yitong Deng
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#include "Points.h"

namespace Meso {
	void Points::Resize(const int _size)
	{
		size = _size;
		for (auto iter = att_map.begin(); iter != att_map.end(); iter++) {
			(iter->second)->Resize(size);
		}
	}

	int Points::Size(void) const
	{
		return size;
	}
}