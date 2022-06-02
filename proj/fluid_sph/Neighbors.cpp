//////////////////////////////////////////////////////////////////////////
// Neighbors
// Copyright (c) (2022-), Yitong Deng
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#include "Neighbors.h"

namespace Meso {

	void Neighbors::set(const Array<int>& _nbs, const Array<int>& _b_nbs, real _r)
	{
		nbs = _nbs;
		b_nbs = _b_nbs;
		r = _r;
	}

}