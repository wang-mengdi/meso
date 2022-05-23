//////////////////////////////////////////////////////////////////////////
// Neighbors
// Copyright (c) (2022-), Bo Zhu, Mengdi Wang, Yitong Deng
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#include "Neighbors.h"

namespace Meso {
	Neighbors::Neighbors()
	{
		e_nbs.clear();
		l_nbs.clear();
	}

	void Neighbors::set_e(const Array<int>& nbs)
	{
		e_nbs = nbs;
	}

	void Neighbors::set_l(const Array<int>& nbs)
	{
		l_nbs = nbs;
	}

	int Neighbors::size_e(void)const
	{
		return e_nbs.size();
	}

	int Neighbors::size_l(void)const
	{
		return l_nbs.size();
	}

	int Neighbors::e(int k)const
	{
		return e_nbs[k];
	}

	int Neighbors::l(int k)const
	{
		return l_nbs[k];
	}
}