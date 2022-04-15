//////////////////////////////////////////////////////////////////////////
// Face data on MacGrid
// Copyright (c) (2022-), Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "Grid.h"

namespace Meso {

	template<class T, int d, DataHolder side = DataHolder::HOST>
	class FaceField {
		Typedef_VectorD(d);
	public:
		Grid<d, CENTER> grid;
		//std::array<Array<T, side>, d> face_data;
		std::shared_ptr<Array<T, side>> face_data[3] = { nullptr,nullptr,nullptr };
		FaceField() {}
		FaceField(const Grid<d, CENTER>& _grid) { Init(_grid); }
		FaceField(const Grid<d, CENTER>& _grid, const T value) { Init(_grid);  Fill(value); }
		void Fill(const T value) { for (int axis = 0; axis < d; axis++) ArrayFunc::Fill(*face_data[axis], value); }
		void Init(const Grid<d, CENTER>& _grid) {
			grid = _grid;
			for (int axis = 0; axis < d; axis++) {
				int n = grid.Face_DoF(axis);
				if (face_data[axis] == nullptr) face_data[axis] = std::make_shared<Array<T, side>>(n);
				else face_data[axis]->resize(n);
				checkCudaErrors(cudaGetLastError());
			}
		}

		template<DataHolder side1> 
		void Deep_Copy(const FaceField<T, d, side1>& f1) {
			for (int i = 0; i < d; i++) { ArrayFunc::Copy(*face_data[i], f1.Data(i)); }
		}

		inline T& operator()(const int axis, const VectorDi face) { return (*(face_data[axis]))[grid.Face_Index(axis, face)]; }
		inline const T& operator()(int axis, const VectorDi face) const { return (*(face_data[axis]))[grid.Face_Index(axis, face)]; }

		constexpr Array<T, side>& Data(const int axis)noexcept { return *face_data[axis]; }
		constexpr const Array<T, side>& Data(const int axis)const noexcept { return *face_data[axis]; }
		constexpr T* Data_Ptr(const int axis) noexcept { return thrust::raw_pointer_cast(face_data[axis]->data()); }
		constexpr const T* Data_Ptr(const int axis) const noexcept { return thrust::raw_pointer_cast(face_data[axis]->data()); }

		template<class IFFunc>
		void Iterate_Faces(IFFunc f) {
			for (int axis = 0; axis < d; axis++) {
				int n = grid.Face_DoF(axis);
				for (int i = 0; i < n; i++) {
					VectorDi face = grid.Face_Coord(axis, i);
					f(axis, face);
				}
			}
		}

		template<class IFFuncT>
		void Calc_Faces(IFFuncT f) {
			for (int axis = 0; axis < d; axis++) {
				const int dof = grid.Face_DoF(axis);
				thrust::counting_iterator<int> idxfirst(0);
				thrust::counting_iterator<int> idxlast = idxfirst + dof;
				thrust::transform(
					idxfirst,
					idxlast,
					face_data[axis]->begin(),
					[f, axis, this](const int idx) {
						return f(axis, grid.Face_Coord(axis, idx));
					}
				);
			}
		}
	};

	template<class T, int d> using FaceFieldDv = FaceField<T, d, DEVICE>;

}