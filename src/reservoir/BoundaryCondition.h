//////////////////////////////////////////////////////////////////////////
// Boundary Conditions
// Copyright (c) (2022-), Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once

#include "Field.h"
#include "FaceField.h"
#include <thrust/iterator/constant_iterator.h>

namespace Meso {
	template<class DataStructure>
	class BoundaryConditionBase {
	public:
		virtual void Copy_Masked(DataStructure& dest, const DataStructure& src) = 0;
	};

	template<class T, int d, DataHolder side>
	class BoundaryConditionDirect : BoundaryConditionBase<Field<T, d, side>> {
	public:
		Field<bool, d, side> fixed;
		template<DataHolder side1>
		void Init(const Field<bool, d, side1>& _fixed) {
			fixed = _fixed;
		}
		void Copy_Masked(Field<T, d, side>& dest, const Field<T, d, side>& src) {
			ArrayFunc::Copy_Masked(dest.Data(), src.Data(), fixed.Data());
		}
		void Set_Masked(Field<T, d, side>& dest, const T val) {
			ArrayFunc::Copy_Masked(dest.Data(), thrust::make_constant_iterator(val), fixed.Data());
		}
	};
}