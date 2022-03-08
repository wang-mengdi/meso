//////////////////////////////////////////////////////////////////////////
// Basic grid data (with data included)
// Copyright (c) (2022-), Bo Zhu, Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "Grid.h"
#include <thrust/iterator/zip_iterator.h>
#include <thrust/iterator/counting_iterator.h>

namespace Meso {
	template<class T, int d, DataHolder side = DataHolder::HOST>
	class Field {
		Typedef_VectorD(d);
	public:
		Grid<d, GridType::CELL> grid;
		Array<T, side> data;
		Field() {}
		Field(const Grid<d, GridType::CELL>& _grid) { Init(_grid); }
		Field(const Grid<d, GridType::CELL>& _grid, const T value) { Init(_grid, value); }
		void Init(const Grid<d, GridType::CELL>& _grid) {
			grid = _grid;
			data.resize(grid.DoF());
		}
		void Init(const Grid<d, GridType::CELL>& _grid, const T value) {
			Init(_grid);
			ArrayFunc::Fill(data, value);
		}
		
		template<DataHolder side1> void Copy(const Field<T, d, side1>& f1) { ArrayFunc::Copy(data, f1.data); }

		inline T& operator()(const VectorDi& coord) { return data[grid.Index(coord)]; }
		inline const T& operator()(const VectorDi& coord) const { return data[grid.Index(coord)]; }

		template<class CFunc>
		void Iterate_Cells(CFunc f){
			const int dof = grid.DoF();
			for (int c = 0; c < dof; c++) {
				f(grid.Coord(c));
			}
		}

		template<class CFuncT>//CFuncT is a function: VectorDi->T, takes the cell index
		void Calc_Cells(CFuncT f) const {
			const int dof = grid.DoF();
			thrust::counting_iterator<int> idxfirst(0);
			thrust::counting_iterator<int> idxlast = idxfirst + dof;
			thrust::transform(
				idxfirst,
				idxlast,
				data.begin(),
				[f, this](const int idx) {
					return f(grid.Coord(idx));
				}
			);
		}
	};

	template<class T, int d> using FieldDv = Field<T, d, DEVICE>;

}