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

		Setup_Attribute(m, real, 1.); //mass
		Setup_Attribute(p, real, 0.); //pressure
		Setup_Attribute(rho, real, 1.);
		Setup_Attribute(nden, real, 1.);
		Setup_Attribute(V, real, 1.); //control volume (or control area in 2D)
		Setup_Attribute(u, VectorD, VectorD::Zero()); //velocity
		Setup_Attribute(acc, VectorD, VectorD::Zero()); //acceleration
		Setup_Attribute(B, int, 0); //whether is boundary

		Setup_Attribute(bnd_phi, real, 0.); //signed distance to boundary
		Setup_Attribute(bnd_n, VectorD, VectorD::Unit(0)); //normal vector to boundary

		Array<Neighbors> nbs_info;

		real dx = 1.;
	};
}