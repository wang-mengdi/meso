//////////////////////////////////////////////////////////////////////////
// Simulator Driver
// Copyright (c) (2022-), Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "Common.h"
#include "MetaData.h"

namespace Meso {
	class Optimizer {
	public:
		virtual void Output(const bf::path base_path, const int iter) = 0;
		virtual void Optimize(OptimizerDriverMetaData& meta_data) = 0;
		virtual bool Is_Converged(OptimizerDriverMetaData& meta_data) = 0;
	};
}