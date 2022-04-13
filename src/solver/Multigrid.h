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

namespace Meso {

	template<class T>
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
		virtual int XDof() const {
			return dof;
		}

		virtual int YDof() const {
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

			//upstroke (coarse->fine)
			for (int i = L - 1; i >= 0; i--) {
				prolongators[i]->Apply(rs[i], xs[i + 1]);
				ArrayFunc::Add(xs[i], rs[i]);
				mappings[i]->Residual(rs[i], xs[i], bs[i]);
				//use bs to temporarily store the data
				postsmoothers[i]->Apply(bs[i], rs[i]);
				ArrayFunc::Add(xs[i], bs[i]);
			}

			ArrayFunc::Copy(x0, xs[0]);
			checkCudaErrors(cudaGetLastError());
		}

		template<int d>
		void Init_Poisson(const PoissonMapping<T, d>& poisson, const int pre_iter = 2, const int post_iter = 2) {
			Check_Cuda_Memory("begin to poisson");

			Typedef_VectorD(d);
			using PoissonPtr = std::shared_ptr<PoissonMapping<T, d>>;

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

			Check_Cuda_Memory("begin to allocate mapping");

			Info("L={}", L);


			////mappings
			mappings.resize(L + 1);
			mappings[0] = std::make_shared<PoissonMapping<T, d>>(poisson);
			for (int i = 1; i <= L; i++) {
				Info("allocate mapping layer {}", i);
				Check_Cuda_Memory("begin to allocate mapping layer");
				mappings[i] = std::make_shared<PoissonMapping<T, d>>(grids[i]);
				PoissonPtr poisson_fine = std::dynamic_pointer_cast<PoissonMapping<T, d>>(mappings[i - 1]);
				PoissonPtr poisson_coarse = std::dynamic_pointer_cast<PoissonMapping<T, d>>(mappings[i]);
				Coarsener<d>::Apply(*poisson_coarse, *poisson_fine);
			}


			Check_Cuda_Memory("begin to allocate restrictors");

			//restrictors
			restrictors.resize(L);
			for (int i = 0; i < L; i++) {
				//i is fine and i+1 is coarse
				restrictors[i] = std::make_shared<Restrictor<T, d>>(grids[i + 1], grids[i]);
			}

			//prolongators
			prolongators.resize(L);
			for (int i = 0; i < L; i++) {
				//i is fine and i+1 is coarse
				prolongators[i] = std::make_shared<Prolongator<T, d>>(grids[i], grids[i + 1]);
			}

			Check_Cuda_Memory("begin to allocate presmoothers");

			//presmoothers
			presmoothers.resize(L);
			for (int i = 0; i < L; i++) {
				PoissonPtr poisson = std::dynamic_pointer_cast<PoissonMapping<T, d>>(mappings[i]);
				presmoothers[i] = std::make_shared<DampedJacobiSmoother<T>>(*poisson, pre_iter, 2.0 / 3.0);
			}

			//postsmoothers
			postsmoothers.resize(L);
			for (int i = 0; i < L; i++) {
				PoissonPtr poisson = std::dynamic_pointer_cast<PoissonMapping<T, d>>(mappings[i]);
				postsmoothers[i] = std::make_shared<DampedJacobiSmoother<T>>(*poisson, post_iter, 2.0 / 3.0);
			}

			Check_Cuda_Memory("begin to allocate direct solver");

			//direct_solver
			DenseMatrixMapping<T> dense_mapping;
			dense_mapping.Init_PoissonLike(grids[L], *mappings[L]);
			direct_solver = std::make_shared<LUDenseSolver<T>>(dense_mapping);
			
			Check_Cuda_Memory("begin to allocate auxillary arrays");

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

			Check_Cuda_Memory("poisson done");
		}
	};

}