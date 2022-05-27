//////////////////////////////////////////////////////////////////////////
// Initializer of a Fluid SPH System
// Copyright (c) (2022-), Bo Zhu, Mengdi Wang, Yitong Deng
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once

#include "FluidSPH.h"
#include "GeometryPrimitives.h"
#include "Json.h"

namespace Meso {
	template<int d>
	class  SPHInitializer{
	public:
		Typedef_VectorD(d); Typedef_MatrixD(d);

		real scale;
		real fineness;

		void Apply(json& j, FluidSPH<d>& fluid) {
			scale = Json::Value(j, "scale", 10.);
			fineness = Json::Value(j, "fineness", 10.);
			int test = Json::Value(j, "test", 0);
			switch (test) {
			case 0: Case_0(j, fluid); break;
			default:Assert(false, "test {} not exist", test); break;
			}
		}

		// Case 0
		void Case_0(json& j, FluidSPH<d>& fluid) {
			std::cout << "Initializing case 0." << std::endl;
			if constexpr (d == 2) {
				fluid.particles.dx = Initialize_Box_Points_2D(Vector2(0,0), Vector2i(10, 10), Vector2(0.1, 0.1), fluid.particles, fluid.particles.xRef());
				std::function<void(const int idx)> instruction = [&](const int idx) {
					fluid.particles.B(idx) = 1;
				};
				fluid.particles.dx = Initialize_Box_Rim_Points_2D(Vector2(0, 0), Vector2i(3, 3), Vector2i(10, 10), Vector2(0.1, 0.1), fluid.particles, fluid.particles.xRef(), true, instruction);
			}
			fluid.Init();
		}
	};
}