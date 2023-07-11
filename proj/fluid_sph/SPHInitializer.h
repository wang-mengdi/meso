//////////////////////////////////////////////////////////////////////////
// Initializer of a Fluid SPH System
// Copyright (c) (2022-), Bo Zhu, Mengdi Wang, Yitong Deng
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once

#include "FluidSPH.h"
#include "ImplicitManifold.h"
#include "AnalyticalBoundary.h"
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
				// init particles
				int s = (int)scale;
				int r = (int)(scale * 0.5);
				int l = s - r;
				fluid.particles.dx = Initialize_Box_Points_2D(Vector2(l,0), Vector2i(r, s), Vector2(1., 1.), fluid.particles, fluid.particles.xRef());
				fluid.b_particles.dx = Initialize_Box_Rim_Points_2D(Vector2(0, 0), Vector2i(3, 3), Vector2i(s, s), Vector2(1., 1.), fluid.b_particles, fluid.b_particles.xRef());
				// init boundary
				fluid.boundary.Add_Obstacle(std::make_shared<Box<d>>(VectorD::Zero(), scale * VectorD::Ones()));
			}
			else if constexpr (d == 3) {
				int s = (int)scale;
				int r = (int)(scale * 0.5);
				int l = s - r;
				fluid.particles.dx = Initialize_Box_Points_3D(Vector3(l, 0, 0), Vector3i(r, s, s), Vector3(1., 1., 1.), fluid.particles, fluid.particles.xRef());
				fluid.b_particles.dx = Initialize_Box_Rim_Points_3D(Vector3(0, 0, 0), Vector3i(3, 3, 3), Vector3i(s, s, s), Vector3(1., 1., 1.), fluid.b_particles, fluid.b_particles.xRef());
				// init boundary
				fluid.boundary.Add_Obstacle(std::make_shared<Box<d>>(VectorD::Zero(), scale * VectorD::Ones()));
			}
			fluid.Init();
		}
	};
}