//////////////////////////////////////////////////////////////////////////
// Initializer of a Fluid Euler System
// Copyright (c) (2022-), Bo Zhu, Mengdi Wang, Yitong Deng
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once

#include "MELP.h"
#include "GeometryPrimitives.h"
#include "Json.h"

namespace Meso {
	template<int d>
	class  MELPInitializer{
	public:
		Typedef_VectorD(d); Typedef_MatrixD(d);

		real scale;
		real fineness;

		void Apply(json& j, MELP<d>& fluid) {
			scale = Json::Value(j, "scale", 10.);
			fineness = Json::Value(j, "fineness", 10.);
			int test = Json::Value(j, "test", 0);
			switch (test) {
			case 0: Case_0(j, fluid); break;
			default:Assert(false, "test {} not exist", test); break;
			}
		}

		// Case 0
		void Case_0(json& j, MELP<d>& fluid) {
			std::cout << "Initializing case 0." << std::endl;
			if constexpr (d == 3) {
				real dx = fineness;
				fluid.dx = dx;
				Initialize_Lattice_Points(Vector3::Zero(), scale, scale, dx * Vector3::Unit(0), dx * Vector3::Unit(2), fluid.e_particles, "x");
			}
		}
	};
}