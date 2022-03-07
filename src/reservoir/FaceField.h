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
		Grid<d, GridType::CELL> grid;
		std::array<Array<T, side>, d> face_data;
		FaceField() {}
		FaceField(const Grid<d, GridType::CELL>& _grid) { Init(_grid); }
		FaceField(const Grid<d, GridType::CELL>& _grid, const T value) { Init(_grid); for (int axis = 0; axis < d; axis++) ArrayFunc::Fill(face_data[axis], value); }
		void Init(const Grid<d, GridType::CELL>& _grid) {
			grid = _grid;
			for (int axis = 0; axis < d; axis++) face_data[axis].resize(grid.Face_DoF(axis));
		}

		inline T& operator()(const int axis, const VectorDi face) { return face_data[axis][grid.Face_Index(axis, face)]; }
		inline const T& operator()(int axis, const VectorDi face) const { return face_data[axis][grid.Face_Index(axis, face)]; }

		template<class IFFunc>
		void Iterate_Faces(IFFunc f) {
			for (int axis = 0; axis < d; axis++) {
				int n = grid.Face_Dof(axis);
				for (int i = 0; i < n; i++) {
					VectorDi face = grid.Face_Coord(axis, i);
					f(axis, face);
				}
			}
		}
	};

	template<class T, int d> using FaceFieldDv = FaceField<T, d, DEVICE>;

}