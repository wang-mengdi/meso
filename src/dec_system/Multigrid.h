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
#include "Restrictor.h"
#include "Prolongator.h"
#include "DampedJacobiSmoother.h"
#include "GridGSSmoother.h"
#include "LUDirectSolver.h"
#include "SparseDirectSolver.h"
#include "Timer.h"

namespace Meso {

	template<class T, class Restrictor, class Prolongator>
	class VCycleMultigrid :public LinearMapping<T> {
		using LinearMappingPtr = std::shared_ptr<LinearMapping<T>>;
	public:
		int L, dof;

		Array<LinearMappingPtr> mappings;
		Array<LinearMappingPtr> restrictors;//restrictor[i] is applied between layer i and i+1
		Array<LinearMappingPtr> prolongators;//prolongator[i] is applied between layer i and i+1
		Array<LinearMappingPtr> presmoothers;
		Array<LinearMappingPtr> postsmoothers;
		LinearMappingPtr direct_solver;

		Array<ArrayDv<T>> xs;
		Array<ArrayDv<T>> bs;
		Array<ArrayDv<T>> rs;
	public:
		virtual int XDoF() const {
			return dof;
		}

		virtual int YDoF() const {
			return dof;
		}

		//solve some Ax=b
		virtual void Apply(ArrayDv<T>& x0, const ArrayDv<T>& b0) {
			//V-cycle
			//downstroke (fine->coarse)
			ArrayFunc::Copy(bs[0], b0);

			for (int i = 0; i < L; i++) {
				presmoothers[i]->Apply(xs[i], bs[i]);
				mappings[i]->Residual(rs[i], xs[i], bs[i]);
				restrictors[i]->Apply(bs[i + 1], rs[i]);
			}

			checkCudaErrors(cudaGetLastError());
			cudaDeviceSynchronize();
			//direct solve
			direct_solver->Apply(xs[L], bs[L]);
			checkCudaErrors(cudaGetLastError());

			//Assert(ArrayFunc::Is_Finite<T, DEVICE>(xs[L]), "Multigrid error: at coarsest level solved xs[{}]={}", L, xs[L]);

			//upstroke (coarse->fine)
			for (int i = L - 1; i >= 0; i--) {
				prolongators[i]->Apply(rs[i], xs[i + 1]);
				ArrayFunc::Add(xs[i], rs[i]);
				mappings[i]->Residual(rs[i], xs[i], bs[i]);
				//Assert(ArrayFunc::Is_Finite<T, DEVICE>(rs[i]), "Multigrid error: at upstroke level {} rs[i]={}", i, rs[i]);
				//use bs to temporarily store the data
				postsmoothers[i]->Apply(bs[i], rs[i]);
				ArrayFunc::Add(xs[i], bs[i]);
			}

			ArrayFunc::Copy(x0, xs[0]);
			checkCudaErrors(cudaGetLastError());
		}

		//Will add epsilon*I to the system of the coarsest level
		//To make a Poisson system truly positive definite
		template<int d>
		void Init_Poisson(const MaskedPoissonMapping<T, d>& poisson, const int pre_iter = 2, const int post_iter = 2, const T coarsest_add_epsilon = 0) {
			Typedef_VectorD(d);
			using PoissonPtr = std::shared_ptr<MaskedPoissonMapping<T, d>>;

			VectorDi grid_size = poisson.Grid().counts;
			int grid_min_size = grid_size.minCoeff();
			L = (int)std::ceil(log2(grid_min_size)) - 3;
			Array<Grid<d>> grids(L + 1);
			grids[0] = poisson.Grid();
			dof = grids[0].DoF();
			for (int i = 1; i <= L; i++) {
				grid_size /= 2;
				grids[i] = Grid<d>(grid_size);
			}

			////mappings
			mappings.resize(L + 1);
			mappings[0] = std::make_shared<MaskedPoissonMapping<T, d>>(poisson);
			for (int i = 1; i <= L; i++) {
				mappings[i] = std::make_shared<MaskedPoissonMapping<T, d>>(grids[i]);
				PoissonPtr poisson_fine = std::dynamic_pointer_cast<MaskedPoissonMapping<T, d>>(mappings[i - 1]);
				PoissonPtr poisson_coarse = std::dynamic_pointer_cast<MaskedPoissonMapping<T, d>>(mappings[i]);
				Coarsener<d>::Apply(*poisson_coarse, *poisson_fine);
			}

			//restrictors
			restrictors.resize(L);
			for (int i = 0; i < L; i++) {
				//i is fine and i+1 is coarse
				restrictors[i] = std::make_shared<Restrictor>(grids[i + 1], grids[i]);
			}

			//prolongators
			prolongators.resize(L);
			for (int i = 0; i < L; i++) {
				//i is fine and i+1 is coarse
				prolongators[i] = std::make_shared<Prolongator>(grids[i], grids[i + 1]);
			}

			//presmoothers and postsmoothers
			presmoothers.resize(L);
			postsmoothers.resize(L);
			for (int i = 0; i < L; i++) {
				PoissonPtr poisson = std::dynamic_pointer_cast<MaskedPoissonMapping<T, d>>(mappings[i]);

				ArrayDv<T> poisson_diag; PoissonLike_Diagonal(poisson_diag, *poisson);
				presmoothers[i] = std::make_shared<DampedJacobiSmoother<T>>(*(mappings[i]), poisson_diag, pre_iter, (T)(2.0 / 3.0));
				postsmoothers[i] = std::make_shared<DampedJacobiSmoother<T>>(*(mappings[i]), poisson_diag, post_iter, (T)(2.0 / 3.0));
				
				//presmoothers[i] = std::make_shared<GridGSSmoother<T, d>>(*poisson, pre_iter, 0);
				//postsmoothers[i] = std::make_shared<GridGSSmoother<T, d>>(*poisson, post_iter, 1);
			}

			//direct_solver

			direct_solver = std::make_shared<CholeskySparseSolver<T>>(SparseMatrix_From_PoissonLike(grids[L], *mappings[L], coarsest_add_epsilon));

			//DenseMatrixMapping<T> dense_mapping;
			//DenseMatrixMapping_From_Poisson_Like(dense_mapping, grids[L], *mappings[L], coarsest_add_epsilon);
			//direct_solver = std::make_shared<LUDenseSolver<T>>(dense_mapping);

			//PoissonPtr last_layer_poisson = std::dynamic_pointer_cast<MaskedPoissonMapping<T, d>>(mappings[L]);
			//direct_solver = std::make_shared<GridGSSmoother<T, d>>(*last_layer_poisson, 5);

			//auxillary arrays
			xs.resize(L + 1);
			bs.resize(L + 1);
			rs.resize(L + 1);//seem rs only need L
			for (int i = 0; i <= L; i++) {
				int n = grids[i].DoF();
				xs[i].resize(n);
				bs[i].resize(n);
				rs[i].resize(n);
			}
		}
	};

	template<class T, int d> using VCycleMultigridIntp = VCycleMultigrid<T, RestrictorIntp<T, d>, ProlongatorIntp<T, d>>;
	template<class T, int d> using VCycleMultigridSum = VCycleMultigrid<T, RestrictorSum<T, d>, ProlongatorSum<T, d>>;
}