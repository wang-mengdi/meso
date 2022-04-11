//////////////////////////////////////////////////////////////////////////
// Coarsening operator for multigrid geometric preconditioned conjugate gradient
// Copyright (c) (2022-), Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "AuxFunc.h"
#include "Field.h"
#include "Interpolation.h"
#include "PoissonMapping.h"

namespace Meso {
	template<int d>
	__global__ void Coarsen_Fixed_Kernel(const Grid<d> grid_coarser, bool* coarser_fixed, const Grid<d> grid_finer, const bool* finer_fixed) {
		Typedef_VectorD(d);
		static const int dx[8] = { 0,1,0,1,0,1,0,1 };
		static const int dy[8] = { 0,0,1,1,0,0,1,1 };
		static const int dz[8] = { 0,0,0,0,1,1,1,1 };
		VectorDi coarser_coord = GPUFunc::Thread_Coord<d>(blockIdx, threadIdx);
		bool fixed = true;//default value of fixed
		for (int s = 0; s < (1 << d); s++) {
			VectorDi finer_coord = coarser_coord * 2 + VectorFunc::Vi<d>(dx[s], dy[s], dz[s]);
			//if (grid_finer.Valid(finer_coord)) fixed &= finer_fixed[grid_finer.Index(finer_coord)];
			fixed &= grid_finer.Valid(finer_coord) ? finer_fixed[grid_finer.Index(finer_coord)] : true;
		}
		coarser_fixed[grid_coarser.Index(coarser_coord)] = fixed;
	}

	template<class T, int d>
	__global__ void Coarsen_Vol_Kernel(const int axis, const Grid<d> coarser_grid, T* coarser_data, const Grid<d> finer_grid, const T* finer_data) {
		Typedef_VectorD(d);
		VectorDi coarser_face = GPUFunc::Thread_Coord<d>(blockIdx, threadIdx);
		VectorDi finer_face = coarser_face * 2;
		VectorD finer_frac = VectorD::Ones() * 0.5;
		finer_frac[axis] = 0;
		coarser_data[coarser_grid.Index(coarser_face)] = Interpolation::Linear_Intp(finer_grid, finer_data, finer_face, finer_frac);
	}

	template<int d>
	class Coarsener {
		Typedef_VectorD(d);
	public:

		//coarser_poisson must be already allocated previously
		template<class T>
		static void Apply(PoissonMapping<T, d>& coarser_poisson, const decltype(coarser_poisson) finer_poisson) {
			const auto& finer_grid = finer_poisson.fixed.grid;
			const auto& coarser_grid = coarser_poisson.fixed.grid;
			
			//fill fixed
			bool* coarser_fixed = coarser_poisson.fixed.Data();
			const bool* finer_fixed = finer_poisson.fixed.Data();
			coarser_grid.Exec_Kernel(&Coarsen_Fixed_Kernel<d>, coarser_grid, coarser_fixed, finer_grid, finer_fixed);

			//fill vol
			for (int i = 0; i < d; i++) {
				
			}
		}
	};
}