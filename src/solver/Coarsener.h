//////////////////////////////////////////////////////////////////////////
// Coarsening operator for multigrid geometric preconditioned conjugate gradient
// Copyright (c) (2022-), Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "Field.h"

namespace Meso {
	template<class T, int d>
	class Coarsener {
	public:
		void Apply(FieldDv<bool, d>& fixed_coarser, const FieldDv<bool, d>& fixed_finer) {

		}
	};
}