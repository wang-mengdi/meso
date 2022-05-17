//////////////////////////////////////////////////////////////////////////
// Basic grid data representation
// Copyright (c) (2022-), Bo Zhu, Mengdi Wang, Yitong Deng
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "NAParticles.h"
#include "Points.h"

namespace Meso {
	template<int d>
	class EulerParticles : public NAParticles<d> {
		Typedef_VectorD(d); Typedef_MatrixD(d);
	public:
		EulerParticles() : NAParticles<d>()  {
			Init_Attribute_rho();
			Init_Attribute_eta();
			Init_Attribute_nden();
			Init_Attribute_u();
			Init_Attribute_E();
		}

		Setup_Attribute(rho, real, 1.0);
		Setup_Attribute(eta, real, 1.0);
		Setup_Attribute(nden, real, 1.0);
		Setup_Attribute(u, VectorD, VectorD::Zero());
		Setup_Attribute(E, MatrixD, MatrixD::Identity());
	};
}