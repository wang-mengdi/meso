//////////////////////////////////////////////////////////////////////////
// Neighbors
// Copyright (c) (2022-), Bo Zhu, Mengdi Wang, Yitong Deng
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once

#include "Common.h"

namespace Meso {
	class Neighbors {
	private:
		Array<int> e_nbs;
		Array<int> l_nbs;
	public:
		Neighbors();
		void set_e(const Array<int>& nbs);
		void set_l(const Array<int>& nbs);
		int size_e(void)const;
		int size_l(void)const;
		int e(int k)const;
		int l(int k)const;
	};
}