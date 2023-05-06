//////////////////////////////////////////////////////////////////////////
// Simulator Driver
// Copyright (c) (2022-), Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once

#include "Common.h"
#include "MetaData.h"

namespace Meso {
	class Simulator {
	public:
		virtual void Output(DriverMetaData& metadata) = 0;
		virtual void Advance(DriverMetaData& metadata) = 0;
		virtual void Load_Frame(DriverMetaData& metadata) {//load metadata.current_frame
			Error("Simulator::Load_Frame() function not implemented");
		}
		//return: dt, running cfl
		//dt can be inf
		virtual real CFL_Time(const real cfl) = 0;
	};
}