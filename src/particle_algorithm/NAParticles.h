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
	template<int d, class NeighborSearcherType>
	class NAParticles : public Points {
		Typedef_VectorD(d);
	public:
		NeighborSearcherType nbs_searcher;

		NAParticles() {	}

		Setup_Attribute(x, VectorD, VectorD::Zero());
		
	//protected: 
		// 1. Attribute<VectorD> _shit = Attribute<VectorD>(VectorD::Zero());
		// 2. Attribute<VectorD> _shit(VectorD::Zero());
	//public:

		void Update_Searcher(void) {
			nbs_searcher.Update_Points(xRef());
		}
	};
}