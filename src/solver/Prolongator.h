//////////////////////////////////////////////////////////////////////////
// Prolongator in Multigrid
// Copyright (c) (2022-), Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "LinearMapping.h"
#include "Grid.h"
#include "AuxFunc.h"
#include "Interpolation.h"

namespace Meso {
	template<class T, int d>
	__global__ void Prolongator_Kernel(const Grid<d> fine_grid, T* fine_data, const Grid<d> coarse_grid, const T* coarse_data) {
		Typedef_VectorD(d);
		VectorDi fine_coord = GPUFunc::Thread_Coord<d>(blockIdx, threadIdx);
		VectorDi coarse_coord;
		VectorD coarse_frac;
		for (int i = 0; i < d; i++) {
			int x = fine_coord[i];
			if (x & 0x1) {
				coarse_coord[i] = (x / 2);
				coarse_frac[i] = 0.25;
			}
			else {
				coarse_coord[i] = (x / 2) - 1;
				coarse_frac[i] = 0.75;
			}
		}
		fine_data[fine_grid.Index(fine_coord)] = Interpolation::Linear_Intp_Padding0(coarse_grid, coarse_data, coarse_coord, coarse_frac);
	}

	template<class T, int d>
	class Prolongator : public LinearMapping<T> {
		Typedef_VectorD(d);
	public:
		Grid<d> fine_grid, coarse_grid;
		void Init(const Grid<d> _fine, const Grid<d> _coarse) {
			fine_grid = _fine;
			coarse_grid = _coarse;
		}
		//number of cols
		virtual int XDof() const {
			return coarse_grid.DoF();
		}

		//number of rows
		virtual int YDof() const {
			return fine_grid.DoF();
		}

		//input p, get Ap
		virtual void Apply(ArrayDv<T>& fine_data, const ArrayDv<T>& coarse_data) {
			Assert(Size_Match(fine_data, coarse_data), "Prolongator error: mismatch sizes");
			T* fine_ptr = ArrayFunc::Data<T, DEVICE>(fine_data);
			const T* coarse_ptr = ArrayFunc::Data<T, DEVICE>(coarse_data);
			fine_grid.Exec_Kernel(&Prolongator_Kernel<T, d>, fine_grid, fine_ptr, coarse_grid, coarse_ptr);
		}
	};
}