//////////////////////////////////////////////////////////////////////////
// Simulator Driver
// Copyright (c) (2022-), Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once

#include "Common.h"

namespace Meso {
	class Simulator {
	public:
		virtual void Output(const bf::path base_path, const bf::path frame_path) = 0;
		virtual void Advance(const int current_frame, const real current_time, const real dt) = 0;
		//can return inf
		virtual real CFL_Time(const real cfl) = 0;
	};
}