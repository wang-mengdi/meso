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
	//template<class DataStructure>
	//class BoundaryConditionBase {
	//public:
	//	virtual void Copy_Masked(DataStructure& dest, const DataStructure& src) = 0;
	//};

	template<template <class> class BC, class DATA>
	concept BoundaryCondition =
		requires(BC<DATA> bc) {
		//Copy_Masked(DATA& dest, const DATA& src);
		bc.Copy_Masked(std::declval<DATA&>(), std::declval<const DATA&>());
	};

	template<class DataStructure>
	class BoundaryConditionDirect {
	public:
	};

	template<class T, int d, DataHolder side>
	class BoundaryConditionDirect<Field<T, d, side>> {
	public:
		Field<bool, d, side> fixed;
		template<DataHolder side1>
		void Init(const Field<bool, d, side1>& _fixed) {
			fixed = _fixed;
		}
		virtual void Copy_Masked(Field<T, d, side>& dest, const Field<T, d, side>& src) {
			ArrayFunc::Copy_Masked(dest.Data(), src.Data(), fixed.Data());
		}
		virtual void Set_Masked(Field<T, d, side>& dest, const T val) {
			thrust::transform_if(
				dest.Data().begin(),//first
				dest.Data().end(),//last
				fixed.Data().begin(),//stencil
				dest.Data().begin(),//result
				[=]__device__(const T dest_val) { return val; },//op
				thrust::identity<bool>()//pred
			);
		}
	};

	template<class T, int d, DataHolder side>
	class BoundaryConditionDirect<FaceField<T, d, side>> {
		FaceField<bool, d, side> fixed;
		template<DataHolder side1>
		void Init(const FaceField<bool, d, side1>& _fixed) {
			fixed = _fixed;
		}
		virtual void Copy_Masked(FaceField<T, d, side>& dest, const FaceField<T, d, side>& src) {
			for (int axis = 0; axis < d; axis++) {
				ArrayFunc::Copy_Masked(dest.Data(axis), src.Data(axis), fixed.Data(axis));
			}
		}
	};
}