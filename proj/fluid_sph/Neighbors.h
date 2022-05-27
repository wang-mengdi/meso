//////////////////////////////////////////////////////////////////////////
// Neighbors
// Copyright (c) (2022-), Yitong Deng
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once

#include "Common.h"

namespace Meso {
	class Neighbors {
	private:
		Array<int> nbs;
	public:
		Neighbors();
		void set(const Array<int>& nbs);
		int size(void)const;
		int operator[](int k)const;
		const Array<int>& operator()(void)const;
	};
}