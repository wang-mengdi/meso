//////////////////////////////////////////////////////////////////////////
// Basic grid data representation
// Copyright (c) (2022-), Bo Zhu, Mengdi Wang, Yitong Deng
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "Points.h"

namespace Meso {
	template<int d>
	class NAParticles : public Points {
		Typedef_VectorD(d);
	public:
		NAParticles() {
			Add_Attribute<VectorD>("x", VectorD::Zero());
		}
	};
}