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
#include "LUDirectSolver.h"
#include "SparseDirectSolver.h"
#include "Timer.h"

namespace Meso {

	//NOTE: you can use Init_Poisson() to create the system;
	//or first Allocate_Poisson(), then Update_Poisson().
	template<class T, int d, class Restrictor, class Prolongator>
	class VCycleMultigrid :public LinearMapping<T> {
		using LinearMappingPtr = std::shared_ptr<LinearMapping<T>>;
		using MaskedPoissonMappingPtr = std::shared_ptr<MaskedPoissonMapping<T, d>>;
		using DampedJacobiSmootherPtr = std::shared_ptr<DampedJacobiSmoother<T, d>>;
	public:
		int L, dof = 0;
		Array<MaskedPoissonMappingPtr> mappings;
		Array<LinearMappingPtr> restrictors;//restrictor[i] is applied between layer i and i+1
		Array<LinearMappingPtr> prolongators;//prolongator[i] is applied between layer i and i+1
		Array<DampedJacobiSmootherPtr> levelsmoothers;
		DampedJacobiSmootherPtr bottomsmoother;

		Array<ArrayDv<T>> xs;
		Array<ArrayDv<T>> bs;
		Array<ArrayDv<T>> rs;
	public:

		virtual int XDoF() const { return dof; }

		virtual int YDoF() const { return dof; }

		//solve some Ax=b
		virtual void Apply(ArrayDv<T>& x0, const ArrayDv<T>& b0) {
			//V-cycle
			//downstroke (fine->coarse)

			auto add = [=]__device__(T & a, T & b) { a += b; };

			cudaMemset(ArrayFunc::Data(x0), 0, sizeof(T)* x0.size());
			levelsmoothers[0]->Boundary_Apply(x0, b0);
			levelsmoothers[0]->Apply(x0, b0);
			mappings[0]->Residual(rs[0], x0, b0);
			restrictors[0]->Apply(bs[1], rs[0]);
			for (int i = 1; i < L; i++) {
				cudaMemset(ArrayFunc::Data(xs[i]), 0, sizeof(T)* xs[i].size());
				levelsmoothers[i]->Boundary_Apply(xs[i], bs[i]);
				levelsmoothers[i]->Apply(xs[i], bs[i]);
				mappings[i]->Residual(rs[i], xs[i], bs[i]);
				restrictors[i]->Apply(bs[i + 1], rs[i]);
			}

			// bottom
			cudaMemset(ArrayFunc::Data(xs[L]), 0, sizeof(T)* xs[L].size());
			bottomsmoother->Apply(xs[L], bs[L]);

			//upstroke (coarse->fine)
			for (int i = L - 1; i >= 1; i--) {
				prolongators[i]->Apply(rs[i], xs[i + 1]);
				GPUFunc::Cwise_Mapping_Wrapper(ArrayFunc::Data(xs[i]), ArrayFunc::Data(rs[i]), add, xs[i].size());
				levelsmoothers[i]->Apply(xs[i], bs[i]);
				levelsmoothers[i]->Boundary_Apply(xs[i], bs[i]);
			}
			prolongators[0]->Apply(rs[0], xs[1]);
			GPUFunc::Cwise_Mapping_Wrapper(ArrayFunc::Data(x0), ArrayFunc::Data(rs[0]), add, x0.size());
			levelsmoothers[0]->Apply(x0, b0);
			levelsmoothers[0]->Boundary_Apply(x0, b0);
		}

		//Update poisson system to an already allocated system
		template<int d>
		void Update_Poisson(const MaskedPoissonMapping<T, d>& poisson) {
			mappings[0]->Init(poisson.cell_type, poisson.vol);
			for (int i = 1; i <= L; i++) {
				MaskedPoissonMappingPtr poisson_fine = mappings[i - 1];
				MaskedPoissonMappingPtr poisson_coarse = mappings[i];
				Coarsener<d>::Apply(*poisson_coarse, *poisson_fine);
				prolongators[i - 1] = std::make_shared<Prolongator>(poisson_fine->cell_type, poisson_coarse->cell_type);
				//i is fine and i+1 is coarse
				restrictors[i - 1] = std::make_shared<Restrictor>(poisson_coarse->cell_type, poisson_fine->cell_type);
			}

			for (int i = 0; i < L; i++)
				mappings[i]->Search_Boundary();

			for (int i = 0; i < L; i++)
				PoissonLike_One_Over_Diagonal(levelsmoothers[i]->one_over_diag, *(mappings[i]));

			//bottomsmoother
			PoissonLike_One_Over_Diagonal(bottomsmoother->one_over_diag, *(mappings[L]));
		}

		//Will add epsilon*I to the system of the coarsest level
		//To make a Poisson system truly positive definite
		template<int d>
		void Init_Poisson(const MaskedPoissonMapping<T, d>& poisson, const int level_iter = 2, const int boundary_iter = 20, const int bottom_iter = 20){
			Typedef_VectorD(d);
			Grid<d> grid = poisson.Grid();
			Assert(grid.Is_Unpadded(), "Multigrid::Allocate_Poisson error: grid {} padding not allowed", grid);

			VectorDi grid_size = grid.Counts();
			int grid_min_size = grid_size.minCoeff();
			L = (int)std::ceil(log2(grid_min_size)) - 3;
			Array<Grid<d>> grids(L + 1);
			grids[0] = grid;
			dof = grids[0].Counts().prod();
			for (int i = 1; i <= L; i++) {
				grid_size = MathFunc::Round_Up_To_Align<d>(grid_size / 2, Grid<d>::Block_Size());
				grids[i] = Grid<d>(grid_size);
			}

			////mappings
			mappings.resize(L + 1);
			for (int i = 0; i <= L; i++) {
				mappings[i] = std::make_shared<MaskedPoissonMapping<T, d>>(grids[i]);
			}

			//restrictors and prolongators
			restrictors.resize(L);
			prolongators.resize(L);


			//presmoothers and postsmoothers
			levelsmoothers.resize(L);

			//auxillary arrays
			xs.resize(L + 1);
			bs.resize(L + 1);
			rs.resize(L + 1);//seem rs only need L
			for (int i = 0; i <= L; i++) {
				int n = grids[i].Counts().prod();
				if (i != 0)
				{
					xs[i].resize(n);
					bs[i].resize(n);
				}
				if (i != L)
					rs[i].resize(n);
			}

			mappings[0]->Init(poisson.cell_type, poisson.vol);
			for (int i = 1; i <= L; i++) {
				MaskedPoissonMappingPtr poisson_fine = mappings[i - 1];
				MaskedPoissonMappingPtr poisson_coarse = mappings[i];
				Coarsener<d>::Apply(*poisson_coarse, *poisson_fine);
				prolongators[i - 1] = std::make_shared<Prolongator>(poisson_fine->cell_type, poisson_coarse->cell_type);
				//i is fine and i+1 is coarse
				restrictors[i - 1] = std::make_shared<Restrictor>(poisson_coarse->cell_type, poisson_fine->cell_type);
			}

			for (int i = 0; i < L; i++)
				mappings[i]->Search_Boundary();
			
			for (int i = 0; i < L; i++) {
				PoissonLike_One_Over_Diagonal(rs[i], *(mappings[i]));
				levelsmoothers[i] = std::make_shared<DampedJacobiSmoother<T, d>>(*(mappings[i]), rs[i], level_iter, boundary_iter, (T)(2.0 / 3.0));
			}

			//bottomsmoother
			PoissonLike_One_Over_Diagonal(xs[L], *(mappings[L]));
			bottomsmoother = std::make_shared<DampedJacobiSmoother<T, d>>(*(mappings[L]), xs[L], bottom_iter, 0, (T)(2.0 / 3.0));
		}
	};

	template<class T, int d> using VCycleMultigridIntp = VCycleMultigrid<T, d,  RestrictorIntp<T, d>, ProlongatorIntp<T, d>>;
	template<class T, int d> using VCycleMultigridSum = VCycleMultigrid<T, d, RestrictorSum<T, d>, ProlongatorSum<T, d>>;
}