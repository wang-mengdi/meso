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
	__global__ void Coarsen_Cell_Type_Kernel(const GridIndexer<d> grid_coarser, unsigned char* coarser_cell_type, const Grid<d> grid_finer, const unsigned char* finer_cell_type) {
		Typedef_VectorD(d);
		static const int dx[8] = { 0,1,0,1,0,1,0,1 };
		static const int dy[8] = { 0,0,1,1,0,0,1,1 };
		static const int dz[8] = { 0,0,0,0,1,1,1,1 };
		VectorDi coarser_coord = GPUFunc::Thread_Coord<d>(blockIdx, threadIdx);
		unsigned cell_type = 2; //lowest priority for SOLID
		for (int s = 0; s < (1 << d); s++) {
			VectorDi finer_coord = coarser_coord * 2 + MathFunc::Vi<d>(dx[s], dy[s], dz[s]);
			if (grid_finer.Valid(finer_coord))
			{
				if (finer_cell_type[grid_finer.Index(finer_coord)] == 1)
					cell_type = 1;
				else if (finer_cell_type[grid_finer.Index(finer_coord)] == 0 && cell_type != 1)
					cell_type = 0;
			}
		}
		coarser_cell_type[grid_coarser.Index(coarser_coord)] = cell_type;
	}

	template<class T, int d>
	__global__ void Coarsen_Vol_Kernel(const int axis, const GridIndexer<d> coarser_grid, T* coarser_data, const GridIndexer<d> finer_grid, const T* finer_data) {
		Typedef_VectorD(d);
		VectorDi coarser_face = GPUFunc::Thread_Coord<d>(blockIdx, threadIdx);
		VectorDi finer_face = coarser_face * 2;
		VectorD finer_frac = VectorD::Ones() * 0.5;
		finer_frac[axis] = 0;
		coarser_data[coarser_grid.Index(coarser_face)] = IntpLinearPadding0::Value(finer_grid, finer_data, finer_face, finer_frac);
	}

	template<int d>
	class Coarsener {
		Typedef_VectorD(d);
	public:

		//coarser_poisson must be already allocated previously
		template<class T>
		static void Apply(MaskedPoissonMapping<T, d>& coarse_poisson, const decltype(coarse_poisson) fine_poisson) {
			const auto& coarse_grid = coarse_poisson.Grid();
			const auto& fine_grid = fine_poisson.Grid();
						
			//fill fixed
			unsigned char* coarse_cell_type = coarse_poisson.cell_type.Data_Ptr();
			const unsigned char* fine_cell_type = fine_poisson.cell_type.Data_Ptr();
			coarse_grid.Exec_Kernel(&Coarsen_Cell_Type_Kernel<d>, coarse_grid, coarse_cell_type, fine_grid, fine_cell_type);

			//fill vol
			for (int axis = 0; axis < d; axis++) {
				Grid<d> coarse_face_grid = coarse_grid.Face_Grid(axis);
				Grid<d> fine_face_grid = fine_grid.Face_Grid(axis);
				T* coarse_vol = coarse_poisson.vol.Data_Ptr(axis);
				const T* fine_vol = fine_poisson.vol.Data_Ptr(axis);
				coarse_face_grid.Exec_Kernel(&Coarsen_Vol_Kernel<T, d>, axis, coarse_face_grid, coarse_vol, fine_face_grid, fine_vol);
			}
			cudaDeviceSynchronize();
			checkCudaErrors(cudaGetLastError());
		}
	};
}