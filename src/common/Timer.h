//////////////////////////////////////////////////////////////////////////
// Timer
// Copyright (c) (2022-), Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include <chrono>
#include <ctime>
#include "Common.h"
#include "Constants.h"

namespace Meso {

	class Timer {
	public:
		std::chrono::time_point<std::chrono::system_clock> total_start;
		std::chrono::time_point<std::chrono::system_clock> lap_start;
		Timer() {
			Reset();
		}
		void Reset(void);
		//total time
		real Total_Time(const real unit = PhysicalUnits::s);
		//lap time, and reset the lap clock
		real Lap_Time(const real unit = PhysicalUnits::s);
	};

}