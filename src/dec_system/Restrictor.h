//////////////////////////////////////////////////////////////////////////
// Restrictor in Multigrid
// Copyright (c) (2022-), Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "LinearMapping.h"
#include "Grid.h"
#include "AuxFunc.h"

namespace Meso {

	template<class T, int d>
	__global__ void Restrictor_Intp_Axis_Kernel(const int axis, const Grid<d> grid, T* intp_data, const T* original_data) {
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
	__global__ void Restrictor_Intp_Coarser_Kernel(const Grid<d> coarser_grid, T* coarser_data, const Grid<d> finer_grid, const T* finer_data) {
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
	class RestrictorIntp : public LinearMapping<T> {
	public:
		Grid<d> coarser_grid, finer_grid;
		ArrayDv<T> intp_data_old;
		ArrayDv<T> intp_data_new;

		RestrictorIntp() {}
		RestrictorIntp(const Grid<d> coarse, const Grid<d> fine) {
			Init(coarse, fine);
		}

		void Init(const Grid<d> _coarser, const Grid<d> _finer) {
			coarser_grid = _coarser, finer_grid = _finer;
			intp_data_old.resize(XDoF());
			intp_data_new.resize(XDoF());
		}

		//number of cols
		virtual int XDoF() const {
			return finer_grid.DoF();
		}

		//number of rows
		virtual int YDoF() const {
			return coarser_grid.DoF();
		}

		//input p, get Ap
		virtual void Apply(ArrayDv<T>& coarser_data, const ArrayDv<T>& finer_data) {
			Memory_Check(coarser_data, finer_data, "Restrictor::Apply error: not enough space");
			T* intp_ptr_old = ArrayFunc::Data<T, DEVICE>(intp_data_old);
			T* intp_ptr_new = ArrayFunc::Data<T, DEVICE>(intp_data_new);
			const T* original_ptr = ArrayFunc::Data<T, DEVICE>(finer_data);
			for (int axis = 0; axis < d; axis++) {
				finer_grid.Exec_Kernel(&Restrictor_Intp_Axis_Kernel<T, d>, axis, finer_grid, intp_ptr_new, axis == 0 ? original_ptr : intp_ptr_old);
				std::swap(intp_ptr_old, intp_ptr_new);
			}
			coarser_grid.Exec_Kernel(&Restrictor_Intp_Coarser_Kernel<T, d>, coarser_grid, ArrayFunc::Data<T, DEVICE>(coarser_data), finer_grid, intp_ptr_old);
			checkCudaErrors(cudaGetLastError());
		}
	};

	template<class T, int d>
	__global__ void Restrictor_Sum_Coarse_Kernel(const Grid<d> coarse_grid, T* coarse_data, const Grid<d> fine_grid, const T* fine_data) {
		Typedef_VectorD(d);
		VectorDi coarse_coord = GPUFunc::Thread_Coord<d>(blockIdx, threadIdx);
		T sum = 0;
		int bc_num = (1 << d);
		for (int s = 0; s < bc_num; s++) {
			int dx = ((s >> 0) & 1), dy = ((s >> 1) & 1), dz = ((s >> 2) & 1);
			VectorDi fine_coord = coarse_coord * 2 + VectorFunc::Vi<d>(dx, dy, dz);
			if (fine_grid.Valid(fine_coord)) sum += fine_data[fine_grid.Index(fine_coord)];
		}
		coarse_data[coarse_grid.Index(coarse_coord)] = sum / bc_num;
	}

	template<class T, int d>
	class RestrictorSum : public LinearMapping<T> {
	public:
		Grid<d> coarse_grid, fine_grid;

		RestrictorSum() {}
		RestrictorSum(const Grid<d> coarse, const Grid<d> fine) {
			Init(coarse, fine);
		}
		void Init(const Grid<d> _coarse, const Grid<d> _fine) {
			coarse_grid = _coarse, fine_grid = _fine;
		}

		//number of cols
		virtual int XDoF() const { return fine_grid.DoF(); }
		//number of rows
		virtual int YDoF() const { return coarse_grid.DoF(); }

		virtual void Apply(ArrayDv<T>& coarse_data, const ArrayDv<T>& fine_data) {
			Memory_Check(coarse_data, fine_data, "RestrictorSum::Apply error: not enough space");
			T* coarse_ptr = ArrayFunc::Data<T, DEVICE>(coarse_data);
			const T* fine_ptr = ArrayFunc::Data<T, DEVICE>(fine_data);
			coarse_grid.Exec_Kernel(&Restrictor_Sum_Coarse_Kernel<T, d>, coarse_grid, coarse_ptr, fine_grid, fine_ptr);
			checkCudaErrors(cudaGetLastError());
		}
	};
}