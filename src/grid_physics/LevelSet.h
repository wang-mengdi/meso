//////////////////////////////////////////////////////////////////////////
// LevelSet
// Copyright (c) (2022-), Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once

#include "Field.h"

namespace Meso {
	template<class T, int d>
	class LevelSet {
	public:
		Field<T, d> phi;
	};
}