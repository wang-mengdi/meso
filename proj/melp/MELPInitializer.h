//////////////////////////////////////////////////////////////////////////
// Initializer of a Fluid Euler System
// Copyright (c) (2022-), Bo Zhu, Mengdi Wang, Yitong Deng
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once

#include "MELP.h"
//#include "GeometryPrimitives.h"
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
			//if constexpr (d == 3) {
			//	fluid.e_particles.dx = Initialize_Sphere_Points_Regular(Vector3::Zero(), 1., 1000, fluid.e_particles, fluid.e_particles.xRef());
			//	//fluid.e_particles.dx = fineness;
			//	//Initialize_Lattice_Points(Vector3::Zero(), 30, 30, Vector3::Unit(0), Vector3::Unit(1), fluid.e_particles, fluid.e_particles.xRef());
			//	//fluid.e_particles.dx = 1;
			//}
			//fluid.Init();
			if constexpr (d == 3) {
				fluid.e_particles.dx = Initialize_Sphere_Points_Regular(Vector3::Zero(), scale, 1000, fluid.e_particles, fluid.e_particles.xRef());
#pragma omp parallel for
				for (int i = 0; i < fluid.e_particles.Size(); i++) {
					fluid.e_particles.x(i)[0] *= 1.25;
				}
			}
			fluid.Init();
			fluid.enclosed_vol = fluid.Compute_Enclosed_Volume();
			fluid.enclosed_amount = fluid.enclosed_vol * (fluid.p_out);
		}
	};
}