//////////////////////////////////////////////////////////////////////////
// Fluid Euler
// Copyright (c) (2022-), Bo Zhu, Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once

#include "Multigrid.h"
#include "ConjugateGradient.h"
#include "Advection.h"

namespace Meso {
	template<int d>
	class FluidEuler : public Simulator {
		Typedef_VectorD(d);
	public:
		FaceFieldDv<real, d> velocity;
	};
}