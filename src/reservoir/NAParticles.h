//////////////////////////////////////////////////////////////////////////
// Basic grid data representation
// Copyright (c) (2022-), Bo Zhu, Mengdi Wang, Yitong Deng
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "Points.h"
#include "NeighborSearcher.h"
#include "Hashtable.h"

namespace Meso {
	template<int d>
	class NAParticles : public Points {
		Typedef_VectorD(d);
	public:

		std::shared_ptr<NeighborSearcher<d>> nbs_searcher;

		std::map<std::string, std::function<VectorD(const int)>> Get;
		std::map<std::string, std::function<Array<VectorD>&()>> Ref;

		NAParticles() {
			Add_Attribute<VectorD>("x", VectorD::Zero());
			Get["x"] = [&](const int idx)->VectorD {
				return this->template Get_Entry<VectorD>("x", idx);
			};
			Ref["x"] = [&]()->Array<VectorD>& {
				return this->template Get_Attribute<VectorD>("x");
			};
			nbs_searcher = std::make_shared<NeighborKDTree<d>>();
			nbs_searcher->Update_Points(Ref["x"]());
		}
	};
}