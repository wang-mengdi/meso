//////////////////////////////////////////////////////////////////////////
// SPH particle representation
// Copyright (c) (2022-), Bo Zhu, Mengdi Wang, Yitong Deng
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "NAParticles.h"
#include "Neighbors.h"

namespace Meso {
	template<int d>
	class SPHParticles : public NAParticles<d> {
		Typedef_VectorD(d); Typedef_MatrixD(d);
	public:
		SPHParticles() : NAParticles<d>()  {}

		Setup_Attribute(rho, real, 1.0);
		Setup_Attribute(nden, real, 1.0);
		Setup_Attribute(V, real, 1.0); //control volume (or control area in 2D)
		Setup_Attribute(u, VectorD, VectorD::Zero()); //velocity
		Setup_Attribute(B, int, 0); //whether is boundary

		Array<Neighbors> nbs_info;

		real dx = 1.;
	};
}