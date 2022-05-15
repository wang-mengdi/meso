//////////////////////////////////////////////////////////////////////////
// Fluid Euler
// Copyright (c) (2022-), Bo Zhu, Mengdi Wang, Yitong Deng
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once

#include "Simulator.h"
#include "Points.h"

namespace Meso {
	template<int d>
	class MELP : public Simulator {
		Typedef_VectorD(d);
	public:
		void Init() {
			std::cout << "Initializing" << std::endl;
		}
		virtual real CFL_Time(const real cfl) {
			return 0.;
		}
		virtual void Output(const bf::path base_path, const int frame) {
			std::cout << "Writing Output" << std::endl;
		}
		virtual void Advance(const int current_frame, const real current_time, const real dt) {
			std::cout << "Advancing" << std::endl;
		}
	};
}