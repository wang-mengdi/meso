//////////////////////////////////////////////////////////////////////////
// Coarsening operator for multigrid geometric preconditioned conjugate gradient
// Copyright (c) (2022-), Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "AuxFunc.h"
#include "Field.h"
#include "PoissonMapping.h"

namespace Meso {
	template<int d>
	__global__ static void Coarsen_Kernel(const Grid<d> grid_coarser, bool* coarser_fixed, const Grid<d> grid_finer, const bool* finer_fixed) {
		Typedef_VectorD(d);
		static const int dx[8] = { 0,1,0,1,0,1,0,1 };
		static const int dy[8] = { 0,0,1,1,0,0,1,1 };
		static const int dz[8] = { 0,0,0,0,1,1,1,1 };
		VectorDi coarser_coord = VectorFunc::Vi<d>(
			blockIdx.x * grid_coarser.block_size + threadIdx.x,
			blockIdx.y * grid_coarser.block_size + threadIdx.y,
			blockIdx.z * grid_coarser.block_size + threadIdx.z
			);
		bool fixed = true;//default value of fixed
		for (int s = 0; s < (1 << d); s++) {
			VectorDi finer_coord = coarser_coord * 2 + VectorFunc::Vi<d>(dx[s], dy[s], dz[s]);
			//if (grid_finer.Valid(finer_coord)) fixed &= finer_fixed[grid_finer.Index(finer_coord)];
			fixed &= grid_finer.Valid(finer_coord) ? finer_fixed[grid_finer.Index(finer_coord)] : true;
		}
		coarser_fixed[grid_coarser.Index(coarser_coord)] = fixed;
	}

	template<int d>
	class Coarsener {
		Typedef_VectorD(d);
	public:

		template<class T>
		static void Apply(PoissonMapping<T, d>& coarser_poisson, const decltype(coarser_poisson) finer_poisson) {
			const auto& finer_grid = finer_poisson.fixed.grid;
			const auto& coarser_grid = coarser_poisson.fixed.grid;
			
			bool* coarser_fixed = coarser_poisson.fixed.Data();
			const bool* finer_fixed = finer_poisson.fixed.Data();
			coarser_grid.Exec_Kernel(&Coarsen_Kernel<d>, coarser_grid, coarser_fixed, finer_grid, finer_fixed);
			//bool* coarser_fixed = thrust::raw_pointer_cast(coarser_poisson.fixed.data());
			//const bool* finer_fixed=thrust::raw_pointer_cast
		}
		//static void Apply(FieldDv<bool, d>& fixed_coarser, const FieldDv<bool, d>& fixed_finer) {
		//	bool* coarser_data = thrust::raw_pointer_cast(fixed_coarser.data.data());
		//	const bool* finer_data = thrust::raw_pointer_cast(fixed_finer.data.data());
		//	fixed_coarser.grid.Exec_Kernel(&Coarsen_Kernel<d>, fixed_coarser.grid, coarser_data, fixed_finer.grid, finer_data);
		//}
	};
}