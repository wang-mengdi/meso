//////////////////////////////////////////////////////////////////////////
// Basic grid data representation
// Copyright (c) (2022-), Bo Zhu, Mengdi Wang, Yitong Deng
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "Points.h"
#include "NeighborSearcher.h"

namespace Meso {
	template<int d>
	class NAParticles : public Points {
		Typedef_VectorD(d);
	public:

		/*std::shared_ptr<NeighborSearcher<d>> nbs_searcher;*/

		std::map<std::string, std::function<VectorD(const int)>> Get;

		NAParticles() {
			Add_Attribute<VectorD>("x", VectorD::Zero());
			Get["x"] = [&](const int idx)->VectorD {
				return this->template Get_Entry<VectorD>("x", idx);
			};
		}
	};
}