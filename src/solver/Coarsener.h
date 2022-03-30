//////////////////////////////////////////////////////////////////////////
// Coarsening operator for multigrid geometric preconditioned conjugate gradient
// Copyright (c) (2022-), Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "AuxFunc.h"
#include "Field.h"

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

		static void Apply(FieldDv<bool, d>& fixed_coarser, const FieldDv<bool, d>& fixed_finer) {
			bool* coarser_data = thrust::raw_pointer_cast(fixed_coarser.data.data());
			const bool* finer_data = thrust::raw_pointer_cast(fixed_finer.data.data());
			if constexpr (d == 2) {
				int Nx = fixed_coarser.grid.counts[0], Ny = fixed_coarser.grid.counts[1];
				Coarsen_Kernel<d> << <dim3(Nx >> 3, Ny >> 3), dim3(8, 8) >> > (fixed_coarser.grid, coarser_data, fixed_finer.grid, finer_data);
			}
			else if constexpr (d == 3) {
				int Nx = fixed_coarser.grid.counts[0], Ny = fixed_coarser.grid.counts[1], Nz = fixed_coarser.grid.counts[2];
				Coarsen_Kernel<d> << <dim3(Nx >> 2, Ny >> 2, Nz >> 2), dim3(4, 4, 4) >> > (fixed_coarser.grid, coarser_data, fixed_finer.grid, finer_data);
			}
		}
	};
}