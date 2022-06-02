//////////////////////////////////////////////////////////////////////////
// Neighbors
// Copyright (c) (2022-), Yitong Deng
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once

#include "Common.h"

namespace Meso {
	struct Neighbors {
		Array<int> nbs;
		Array<int> b_nbs;
		real r;
		void set(const Array<int>& _nbs, const Array<int>& _b_nbs, real _r);
	};
}