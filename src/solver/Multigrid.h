//////////////////////////////////////////////////////////////////////////
// Multigrid method
// Copyright (c) (2022-), Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "LinearMapping.h"
#include "PoissonMapping.h"
#include "AuxFunc.h"
#include "Coarsener.h"

namespace Meso {

	template<class T, int d>
	class VCycleMultigrid :LinearMapping<T> {
		Typedef_VectorD(d);
	public:
		int L;

		Array<shared_ptr<LinearMapping<T>>> mappings;
		Array<shared_ptr<LinearMapping<T>>> restrictors;
		Array<shared_ptr<LinearMapping<T>>> prolongators;
		Array<shared_ptr<LinearMapping<T>>> presmoothers;
		Array<shared_ptr<LinearMapping<T>>> postsmoothers;
		shared_ptr<LinearMapping<T>> direct_solver;

		ArrayDv<ArrayDv<T>> xs;
		ArrayDv<ArrayDv<T>> bs;
		ArrayDv<ArrayDv<T>> rs;
		ArrayDv<T> x_temp;//store the value calcualted by postsmoother
	public:
		virtual int XDof() const {

		}

		virtual int YDof() const {

		}

		//solve some Ax=b
		virtual void Apply(ArrayDv<T>& x, const ArrayDv<T>& b) {

		}

		void Init_Poisson(const PoissonMapping<T, d> &poisson) {
			VectorDi grid_size = poisson.fixed.grid.counts;
			int grid_min_size = grid_size.minCoeff();
			L = (int)std::ceil(log2(grid_min_size)) - 3;

			//mappings
			mappings.resize(L + 1);
			mappings[0] = std::make_shared<PoissonMapping<T, d>>(poisson);
			for (int i = 1; i <= L; i++) {
				grid_size /= 2;
				mappings[i] = std::make_shared<PoissonMapping<T, d>>();
				mappings[i]->Init(Grid<d>(grid_size));
				Coarsener<d>::Apply(*mappings[i], *mappings[i - 1]);
			}

			//restrictors

		}

		void V_Cycle(ArrayDv<T>& x0, const ArrayDv<T>& b0) {
			ArrayFunc::Copy(b[0], b0);
			for (int i = 0; i < L; i++) {
				presmoothers[i]->Apply(xs[i], b[i]);
				mappings[i]->Residual(rs[i], xs[i], b[i]);
				restrictors[i]->Apply(b[i + 1], rs[i]);
			}
			direct_solver->Apply(xs[L], b[L]);
			for (int i = L - 1; i >= 0; i--) {
				prolongators[i]->Apply(rs[i], xs[i + 1]);
				ArrayFunc::Add(xs[i], rs[i]);
				mappings[i]->Residual(rs[i], xs[i], b[i]);
				postsmoothers[i]->Apply(x_temp, rs[i]);
				ArrayFunc::Add(xs[i], x_temp);
			}
			ArrayFunc::Copy(x0, xs[0]);
		}
	};

}