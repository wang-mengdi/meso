//////////////////////////////////////////////////////////////////////////
// Poisson Linear Mapping
// Copyright (c) (2018-), Zangyueyang Xian, Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "LinearMapping.h"
#include "Field.h"
#include "LambdaHelper.h"

namespace Meso {

	template<class T, int d>
	class PoissonMapping : public LinearMapping<T> {
		Typedef_VectorD(d);
	public:
		int dof;
		FieldDv<T, d> vol;
		FieldDv<bool, d> fixed;

		void Init(const Grid<d, GridType::CELL>& grid, IFFunc<T, d> vol_func, CFunc<T, d> is_unknown_func) {
			dof = grid.DoF();
			vol.Calc_Each(vol_func);
			fixed.Calc_Each(
				[=](const VectorDi& cell)->bool {return !is_unknown_func(cell); }
			);
		}

		virtual int xDoF() const { return dof; }//number of cols

		virtual int yDoF() const { return dof; }//number of rows

		//input p, get Ap
		virtual void applyMapping(ArrayDv<T>& Ap, const ArrayDv<T>& p) {

		}
	};

}