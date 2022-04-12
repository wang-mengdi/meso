//////////////////////////////////////////////////////////////////////////
// Restrictor
// Copyright (c) (2022-), Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "LinearMapping.h"
#include "Grid.h"
#include "AuxFunc.h"

namespace Meso {

	template<class T, int d>
	__global__ void Restrictor_Axis_Kernel(const int axis, const Grid<d> grid, T* intp_data, const T* original_data) {
		Typedef_VectorD(d);
		VectorDi coord = GPUFunc::Thread_Coord<d>(blockIdx, threadIdx);
		T data0, data1, data2, data3;
		int idx = grid.Index(coord);
		data1 = original_data[idx];
		coord[axis]--;
		if (coord[axis] < 0) data0 = 0;
		else data0 = original_data[grid.Index(coord)];
		coord[axis] += 2;
		if (coord[axis] >= grid.counts[axis]) data2 = 0;
		else data2 = original_data[grid.Index(coord)];
		coord[axis]++;
		if (coord[axis] >= grid.counts[axis]) data3 = 0;
		else data3 = original_data[grid.Index(coord)];
		intp_data[idx] = (data0 + 3 * data1 + 3 * data2 + data3) / 8;
	}

	template<class T, int d>
	__global__ void Restrictor_Coarser_Kernel(const Grid<d> coarser_grid, T* coarser_data, const Grid<d> finer_grid, const T* finer_data) {
		Typedef_VectorD(d);
		VectorDi coarser_coord = GPUFunc::Thread_Coord<d>(blockIdx, threadIdx);
		VectorDi finer_coord = coarser_coord * 2;
		int coarser_idx = coarser_grid.Index(coarser_coord);
		if (finer_grid.Valid(finer_coord)) {
			coarser_data[coarser_idx] = finer_data[finer_grid.Index(finer_coord)];
		}
		else coarser_data[coarser_idx] = 0;
	}

	template<class T, int d>
	class Restrictor : public LinearMapping<T> {
	public:
		Grid<d> coarser_grid, finer_grid;
		ArrayDv<T> intp_data;

		void Init(const Grid<d> _coarser, const Grid<d> _finer) {
			coarser_grid = _coarser, finer_grid = _finer;
			intp_data.resize(XDof());
		}

		//number of cols
		virtual int XDof() const {
			return finer_grid.DoF();
		}

		//number of rows
		virtual int YDof() const {
			return coarser_grid.DoF();
		}

		//input p, get Ap
		virtual void Apply(ArrayDv<T>& Ap, const ArrayDv<T>& p) {
			T* intp_ptr = ArrayFunc::Data<T, DEVICE>(intp_data);
			const T* original_ptr = ArrayFunc::Data<T, DEVICE>(p);
			for (int axis = 0; axis < d; axis++) {
				finer_grid.Exec_Kernel(&Restrictor_Axis_Kernel<T, d>, axis, finer_grid, intp_ptr, original_ptr);
			}
			coarser_grid.Exec_Kernel(&Restrictor_Coarser_Kernel<T, d>, coarser_grid, ArrayFunc::Data<T, DEVICE>(Ap), finer_grid, intp_ptr);
		}
	};
}