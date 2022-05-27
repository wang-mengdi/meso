//////////////////////////////////////////////////////////////////////////
// Neighbors
// Copyright (c) (2022-), Yitong Deng
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#include "Neighbors.h"

namespace Meso {
	Neighbors::Neighbors()
	{
		nbs.clear();
	}

	void Neighbors::set(const Array<int>& _nbs)
	{
		nbs = _nbs;
	}

	int Neighbors::size(void)const
	{
		return nbs.size();
	}

	int Neighbors::operator[](int k)const
	{
		return nbs[k];
	}

	const Array<int>& Neighbors::operator()(void)const
	{
		return nbs;
	}

}