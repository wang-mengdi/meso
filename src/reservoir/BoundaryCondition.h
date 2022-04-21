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
	class BoundaryCondition {
	public:
		virtual void Apply(DataStructure& data)const = 0;
		//virtual void Copy_Masked(DataStructure& dest, const DataStructure& src) const = 0;
		//virtual void Copy_UnMasked(DataStructure& dest, const DataStructure& src) const = 0;
	};

	//template<template <class> class BC, class DATA>
	//concept BoundaryCondition =
	//	requires(BC<DATA> bc) {
	//	//Copy_Masked(DATA& dest, const DATA& src);
	//	bc.Copy_Masked(std::declval<DATA&>(), std::declval<const DATA&>());
	//};

	template<class DataStructure>
	class BoundaryConditionDirect : public BoundaryCondition<DataStructure> {
	public:
	};

	template<class T, int d, DataHolder side>
	class BoundaryConditionDirect<Field<T, d, side>> : public BoundaryCondition<Field<T, d, side>> {
	public:
		Array<int, side> indices;
		Array<T, side> values;
		virtual void Apply(Field<T, d, side>& data)const {
			thrust::scatter(
				values.begin(),//first
				values.end(),//last
				indicies.begin(),//map
				data.Data().begin(),//result
				);
		}
		//Field<bool, d, side> fixed;
		//template<DataHolder side1>
		//void Init(const Field<bool, d, side1>& _fixed) {
		//	fixed = _fixed;
		//}
		//virtual void Copy_Masked(Field<T, d, side>& dest, const Field<T, d, side>& src) const {
		//	ArrayFunc::Copy_Masked(dest.Data(), src.Data(), fixed.Data());
		//}
		//virtual void Set_Masked(Field<T, d, side>& dest, const T val) const {
		//	thrust::transform_if(
		//		dest.Data().begin(),//first
		//		dest.Data().end(),//last
		//		fixed.Data().begin(),//stencil
		//		dest.Data().begin(),//result
		//		[=]__device__(const T dest_val) { return val; },//op
		//		thrust::identity<bool>()//pred
		//	);
		//}
	};

	template<class T, int d, DataHolder side>
	class BoundaryConditionDirect<FaceField<T, d, side>> : public BoundaryCondition<FaceField<T, d, side>> {
	public:
		std::array<BoundaryConditionDirect<Field<T, d, side>>, d> field_bc;
		virtual void Apply(FaceField<T, d, side>& face_field)const {
			for (int axis = 0; axis < d; axis++) {
				Field<T, d, side> field(face_field.grid.Face_Grid(axis), face_field.face_data[axis]);
				field_bc.Apply(field);
			}
		}
		//template<DataHolder side1>
		//void Init(const FaceField<bool, d, side1>& _fixed) {
		//	fixed = _fixed;
		//}
		//virtual void Copy_Masked(FaceField<T, d, side>& dest, const FaceField<T, d, side>& src) const {
		//	for (int axis = 0; axis < d; axis++) {
		//		ArrayFunc::Copy_Masked(dest.Data(axis), src.Data(axis), fixed.Data(axis));
		//	}
		//}
		//virtual void Copy_Unmasked(FaceField<T, d, side>& dest, const FaceField<T, d, side>& src)const {
		//	for (int axis = 0; axis < d; axis++) {
		//		ArrayFunc::Copy_UnMasked(dest.Data(axis), src.Data(axis), fixed.Data(axis));
		//	}
		//}
	};
}