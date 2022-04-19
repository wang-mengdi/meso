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
		virtual void Output(const std::string base_path, const std::string frame_path) = 0;
		virtual void Advance(const int current_frame, const real current_time, const real dt) = 0;
	};
}