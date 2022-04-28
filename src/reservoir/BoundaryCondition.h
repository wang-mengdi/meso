//////////////////////////////////////////////////////////////////////////
// Boundary Conditions
// Copyright (c) (2022-), Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once

#include "Field.h"
#include "FaceField.h"
#include <thrust/iterator/constant_iterator.h>
#include <thrust/gather.h>

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
		Typedef_VectorD(d);
	public:
		Array<int, side> indices;
		Array<T, side> values;
		void Init(Field<bool, d>& fixed, Field<T, d> val_field) {
			Array<std::pair<int, T>> bc_pairs;
			fixed.Iterate_Cells(
				[&](const VectorDi cell) {
					if (fixed(cell)) {
						bc_pairs.push_back(std::make_pair(fixed.grid.Index(cell), val_field(cell)));
					}
				}
			);
			std::sort(bc_pairs.begin(), bc_pairs.end());
			int n = bc_pairs.size();
			Array<int> ind_host(n);
			Array<T> val_host(n);
			for (int i = 0; i < n; i++) {
				ind_host[i] = bc_pairs[i].first;
				val_host[i] = bc_pairs[i].second;
			}
			indices = ind_host;
			values = val_host;
		}
		virtual void Apply(Field<T, d, side>& data)const {
			thrust::scatter(
				values.begin(),//first
				values.end(),//last
				indices.begin(),//map
				data.Data().begin()//result
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
		void Init(FaceField<bool, d>& fixed, FaceField<T, d>& value) {
			for (int axis = 0; axis < d; axis++) {
				Grid<d> face_grid = fixed.grid.Face_Grid(axis);
				Field<bool, d> face_fixed(face_grid, fixed.face_data[axis]);
				Field<T, d> face_value(face_grid, value.face_data[axis]);
				field_bc[axis].Init(face_fixed, face_value);
			}
		}
		virtual void Apply(FaceField<T, d, side>& face_field)const {
			for (int axis = 0; axis < d; axis++) {
				Field<T, d, side> field(face_field.grid.Face_Grid(axis), face_field.face_data[axis]);
				field_bc[axis].Apply(field);
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

	template<class DataStructure>
	class BoundaryConditionRefrLinear : public BoundaryCondition<DataStructure> {
	public:
	};

	template<class T, int d, DataHolder side>
	class BoundaryConditionRefrLinear<Field<T, d, side>> : public BoundaryCondition<Field<T, d, side>> {
	public:
		//dst=a*src+b
		Array<int, side> src_list;
		Array<int, side> dst_list;
		Array<T, side> a;
		Array<T, side> b;
		Array<T, side> gathered_data;
		virtual void Apply(Field<T, d, side>& F) {
			thrust::gather(
				src_list.begin(),
				src_list.end(),
				F.Data().begin(),
				gathered_data.begin()
			);
			ArrayFunc::Multiply_Scalar(gathered_data, a);
			ArrayFunc::Add_Scalar(gathered_data, a);
			
		}
	};
}