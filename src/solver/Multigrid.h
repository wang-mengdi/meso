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
	public:
		int L, dof = 0;
		Array<MaskedPoissonMappingPtr> mappings;
		Array<LinearMappingPtr> restrictors;//restrictor[i] is applied between layer i and i+1
		Array<LinearMappingPtr> prolongators;//prolongator[i] is applied between layer i and i+1
		Array<LinearMappingPtr> presmoothers;
		Array<LinearMappingPtr> postsmoothers;
		LinearMappingPtr bottomsmoother;

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
			{
				T* bs0_ptr = thrust::raw_pointer_cast(bs[0].data());
				const T* b0_ptr = thrust::raw_pointer_cast(b0.data());
				cudaMemcpy(bs0_ptr, b0_ptr, sizeof(T) * b0.size(), cudaMemcpyDeviceToDevice);
			}

			auto add = [=]__device__(T & a, T & b) { a += b; };

			for (int i = 0; i < L; i++) {
				cudaMemset(ArrayFunc::Data(xs[i]), 0, sizeof(T)* xs[i].size());
				presmoothers[i]->Apply(xs[i], bs[i]);
				mappings[i]->Residual(rs[i], xs[i], bs[i]);
				restrictors[i]->Apply(bs[i + 1], rs[i]);
			}

			// bottom
			checkCudaErrors(cudaGetLastError());
			cudaDeviceSynchronize();
			bottomsmoother->Apply(xs[L], bs[L]);
			checkCudaErrors(cudaGetLastError());

			//upstroke (coarse->fine)
			for (int i = L - 1; i >= 0; i--) {
				prolongators[i]->Apply(rs[i], xs[i + 1]);
				GPUFunc::Cwise_Mapping_Wrapper(ArrayFunc::Data(xs[i]), ArrayFunc::Data(rs[i]), add, xs[i].size());
				postsmoothers[i]->Apply(xs[i], bs[i]);
			}

			cudaMemcpy(ArrayFunc::Data(x0), ArrayFunc::Data(xs[0]), sizeof(T)* x0.size(), cudaMemcpyDeviceToDevice);

			checkCudaErrors(cudaGetLastError());
		}

		//Allocate most part of the system
		template<int d>
		void Allocate_Poisson(const Grid<d> _grid) {
			Typedef_VectorD(d);
			Assert(_grid.Is_Unpadded(), "Multigrid::Allocate_Poisson error: _grid {} padding not allowed", _grid);

			VectorDi grid_size = _grid.Counts();
			int grid_min_size = grid_size.minCoeff();
			L = (int)std::ceil(log2(grid_min_size)) - 3;
			Array<Grid<d>> grids(L + 1);
			grids[0] = _grid;
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
			presmoothers.resize(L);
			postsmoothers.resize(L);

			//auxillary arrays
			xs.resize(L + 1);
			bs.resize(L + 1);
			rs.resize(L + 1);//seem rs only need L
			for (int i = 0; i <= L; i++) {
				int n = grids[i].Counts().prod();
				xs[i].resize(n);
				bs[i].resize(n);
				rs[i].resize(n);
			}
		}

		template<int d>
		void Update_Poisson_Coarse_Layers(const int level_iter = 2, const int boundary_iter = 20, const int bottom_iter = 20) {

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

			//presmoothers and postsmoothers
			for (int i = 0; i < L; i++) {
				MaskedPoissonMappingPtr poisson = mappings[i];
				ArrayDv<T> poisson_diag; PoissonLike_Diagonal(poisson_diag, *poisson);
				presmoothers[i] = std::make_shared<DampedJacobiSmoother<T, d>>(*(mappings[i]), poisson_diag, level_iter, (T)(2.0 / 3.0));
				postsmoothers[i] = std::make_shared<DampedJacobiSmoother<T, d>>(*(mappings[i]), poisson_diag, level_iter, (T)(2.0 / 3.0));
			}

			//bottomsmoother
			MaskedPoissonMappingPtr poisson = mappings[L];
			ArrayDv<T> poisson_diag; PoissonLike_Diagonal(poisson_diag, *poisson);
			bottomsmoother = std::make_shared<DampedJacobiSmoother<T, d>>(*(mappings[L]), poisson_diag, bottom_iter, (T)(2.0 / 3.0));
		}

		//Update poisson system to an already allocated system
		template<int d>
		void Update_Poisson(const MaskedPoissonMapping<T, d>& poisson, const int level_iter = 2, const int boundary_iter = 20, const int bottom_iter = 20) {
			mappings[0] = std::make_shared<MaskedPoissonMapping<T, d>>(poisson);
			Update_Poisson_Coarse_Layers<d>(level_iter, boundary_iter, bottom_iter);
		}

		//Will add epsilon*I to the system of the coarsest level
		//To make a Poisson system truly positive definite
		template<int d>
		void Init_Poisson(const MaskedPoissonMapping<T, d>& poisson, const int level_iter = 2, const int boundary_iter = 20, const int bottom_iter = 20){
			Allocate_Poisson(poisson.Grid());
			Update_Poisson<d>(poisson, level_iter, boundary_iter, bottom_iter);
		}
	};

	template<class T, int d> using VCycleMultigridIntp = VCycleMultigrid<T, d,  RestrictorIntp<T, d>, ProlongatorIntp<T, d>>;
	template<class T, int d> using VCycleMultigridSum = VCycleMultigrid<T, d, RestrictorSum<T, d>, ProlongatorSum<T, d>>;
}