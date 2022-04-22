//////////////////////////////////////////////////////////////////////////
// Initializer of a Fluid Euler System
// Copyright (c) (2022-), Bo Zhu, Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once

#include "FluidEuler.h"
#include "Json.h"

namespace Meso {
	template<int d>
	class FluidEulerInitializer {
	public:
		void Apply(json& j, FluidEuler<d>& fluid) {
			int test = j.Value("test", 0);
		}
		//an empty field with inflow boundary condition v=1
		void Case_0(json& j, FluidEuler<d>& fluid) {
			
		}
	};
}