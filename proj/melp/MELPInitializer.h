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
		Typedef_VectorD(d);

		void Apply(json& j, MELP<d>& fluid) {
			int test = Json::Value(j, "test", 0);
			switch (test) {
			case 0:Case_0(j, fluid); break;
			default:Assert(false, "test {} not exist", test); break;
			}
		}

		//0
		void Case_0(json& j, MELP<d>& fluid) {
			std::cout << "Initializing case 0." << std::endl;
		}
	};
}