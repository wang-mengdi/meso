//////////////////////////////////////////////////////////////////////////
// Random number
// Copyright (c) (2022-), Bo Zhu, Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once

#include <random>
#include <chrono>
#include "Common.h"

namespace Meso {

	//We will imitate interfaces of Python's random lib here
	namespace Random {
		//int-valued functions
		int RandInt(int a, int b);//random number in [a,b]

		//real-valued functions
		real Random(void);//random number in [0,1)
		real Uniform(real a, real b);
	}

}