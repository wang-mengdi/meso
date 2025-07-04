//////////////////////////////////////////////////////////////////////////
// Simulator Driver
// Copyright (c) (2022-), Fan Feng
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "Common.h"
#include "MetaData.h"

namespace Meso {
	class Optimizer {
	public:
		virtual void Output(OptimizerDriverMetaData& meta_data) = 0;
		virtual void Optimize(OptimizerDriverMetaData& meta_data) = 0;
		virtual bool Is_Converged(OptimizerDriverMetaData& meta_data) = 0;
	};
}