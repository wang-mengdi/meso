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
	__global__ void Coarsen_Vol_Kernel(const int axis, const Grid<d, CORNER> coarser_grid, T* coarser_data, const Grid<d, CORNER> finer_grid, const T* finer_data) {
		Typedef_VectorD(d);
		VectorDi coarser_face = GPUFunc::Thread_Coord<d>(blockIdx, threadIdx);
		VectorDi finer_face = coarser_face * 2;
		VectorD finer_frac = VectorD::Ones() * 0.5;
		finer_frac[axis] = 0;
		coarser_data[coarser_grid.Index(coarser_face)] = Interpolation::Linear_Intp_Padding0(finer_grid, finer_data, finer_face, finer_frac);
	}

	template<int d>
	class Coarsener {
		Typedef_VectorD(d);
	public:

		//coarser_poisson must be already allocated previously
		template<class T>
		static void Apply(PoissonMapping<T, d>& coarse_poisson, const decltype(coarse_poisson) fine_poisson) {
			const auto& coarse_grid = coarse_poisson.Grid();
			const auto& fine_grid = fine_poisson.Grid();
						
			//fill fixed
			bool* coarse_fixed = coarse_poisson.fixed.Data();
			const bool* fine_fixed = fine_poisson.fixed.Data();
			coarse_grid.Exec_Kernel(&Coarsen_Fixed_Kernel<d>, coarse_grid, coarse_fixed, fine_grid, fine_fixed);

			//fill vol
			for (int axis = 0; axis < d; axis++) {
				Grid<d, CORNER> coarse_face_grid = coarse_grid.Face_Grid(axis);
				Grid<d, CORNER> fine_face_grid = fine_grid.Face_Grid(axis);
				T* coarse_vol = coarse_poisson.vol.Data(axis);
				const T* fine_vol = fine_poisson.vol.Data(axis);
				coarse_face_grid.Exec_Kernel(&Coarsen_Vol_Kernel<T, d>, axis, coarse_face_grid, coarse_vol, fine_face_grid, fine_vol);
			}
			cudaDeviceSynchronize();
			checkCudaErrors(cudaGetLastError());
		}
	};
}