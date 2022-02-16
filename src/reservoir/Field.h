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
	template<int d, class CFuncT>
	class Thrust_Unary_Cell_Functor
	{
		Typedef_VectorD(d);
	public:
		using T = decltype(CFuncT(VectorD::Zero()));

		const Grid<d, GridType::CELL>* grid;

		Thrust_Unary_Cell_Functor(const Grid<d, GridType::CELL>* _grid) : grid(_grid) {}

		__host__ __device__
		T operator()(const thrust::tuple<T, int>& x) const {
			int idx = thrust::get<1>(x);
			return CFuncT(grid->Coord(idx));
		}
	};

	template<class T, int d, DataHolder side = DataHolder::HOST>
	class Field {
		Typedef_VectorD(d);
	public:
		Grid<d, GridType::CELL> grid;
		Array<T, side> data;
		Field() {}
		Field(const Grid<d, GridType::CELL>& _grid, const real val = 0) :
			grid(_grid)
		{
			data.resize(grid.DoF());
		}
		inline T& operator()(const VectorDi& coord) { return data[grid.Index(coord)]; }
		inline const T& operator()(const VectorDi& coord) const { return data[grid.Index(coord)]; }

		template<class CFuncT>//CFuncT is a function: VectorDi->T, takes the cell index
		void Calc_Each(CFuncT f) const {
			const int dof = DoF();
			if constexpr (side == DataHolder::HOST) {
#pragma omp parallel for
				for (int c = 0; c < dof; c++) {
					const VectorDi cell = Coord(c);
					data[c] = f(cell);
				}
			}
			else {
				typedef thrust::device_vector<int>::iterator intiter;
				typedef thrust::counting_iterator<int>     countiter;
				typedef thrust::tuple<intiter, countiter>  tpl2intiter;
				typedef thrust::zip_iterator<tpl2intiter>  idxzip;

				thrust::counting_iterator<int> idxfirst(0);
				thrust::counting_iterator<int> idxlast = idxfirst + dof;
				idxzip first = thrust::make_zip_iterator(thrust::make_tuple(data.begin(), idxfirst));
				idxzip  last = thrust::make_zip_iterator(thrust::make_tuple(data.end(), idxlast));
				Thrust_Unary_Cell_Functor<d, CFuncT> my_unary_op;
				thrust::transform(first, last, data.begin(), my_unary_op);
			}
		}
	};

	template<class T, int d> using FieldDv = Field<T, d, DataHolder::DEVICE>;

}