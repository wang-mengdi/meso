//////////////////////////////////////////////////////////////////////////
// Prolongator in Multigrid
// Copyright (c) (2022-), Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "LinearMapping.h"
#include "Grid.h"
#include "AuxFunc.h"

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

			}
			else {
				//coarse_coord[i]=
			}
		}
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

		}
	};
}