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
	__global__ void Restrictor_Intp_Axis_Kernel(const int axis, const GridIndexer<d> grid, T* intp_data, const T* original_data) {
		Typedef_VectorD(d);
		VectorDi coord = GPUFunc::Thread_Coord<d>(blockIdx, threadIdx);
		T data0, data1, data2, data3;
		int idx = grid.Index(coord);
		data1 = original_data[idx];
		coord[axis]--;
		if (coord[axis] < 0) data0 = 0;
		else data0 = original_data[grid.Index(coord)];
		coord[axis] += 2;
		if (coord[axis] >= grid.Dimension(axis)) data2 = 0;
		else data2 = original_data[grid.Index(coord)];
		coord[axis]++;
		if (coord[axis] >= grid.Dimension(axis)) data3 = 0;
		else data3 = original_data[grid.Index(coord)];
		intp_data[idx] = (data0 + 3 * data1 + 3 * data2 + data3) / 8;
	}

	template<class T, int d>
	__global__ void Restrictor_Intp_Coarser_Kernel(const GridIndexer<d> coarser_grid, T* coarser_data, const GridIndexer<d> finer_grid, const T* finer_data) {
		Typedef_VectorD(d);
		VectorDi coarser_coord = GPUFunc::Thread_Coord<d>(blockIdx, threadIdx);
		VectorDi finer_coord = coarser_coord * 2;
		int coarser_idx = coarser_grid.Index(coarser_coord);
		if (finer_grid.Valid(finer_coord)) {
			coarser_data[coarser_idx] = finer_data[finer_grid.Index(finer_coord)] * 4;
		}
		else coarser_data[coarser_idx] = 0;
	}

	template<class T, int d>
	class RestrictorIntp : public LinearMapping<T> {
	public:
		using Base=LinearMapping<T>;
		const FieldDv<unsigned char, d>* coarse_cell_type;
		const FieldDv<unsigned char, d>* fine_cell_type;
		ArrayDv<T> intp_data_old;
		ArrayDv<T> intp_data_new;

		RestrictorIntp() {}
		RestrictorIntp(const FieldDv<unsigned char, d>& coarse_cell_type, const FieldDv<unsigned char, d>& fine_cell_type) {
			Init(coarse_cell_type, fine_cell_type);
		}

		void Init(const FieldDv<unsigned char, d>& _coarse_cell_type, const FieldDv<unsigned char, d>& _fine_cell_type) {
			Assert(_coarse_cell_type.grid.Is_Unpadded(), "RestrictorIntp: _coarser {} invalid, must be unpadded", _coarse_cell_type.grid);
			Assert(_fine_cell_type.grid.Is_Unpadded(), "RestrictorIntp: _finer {} invalid, must be unpadded", _fine_cell_type.grid);
			coarse_cell_type = &_coarse_cell_type, fine_cell_type = &_fine_cell_type;
			intp_data_old.resize(XDoF());
			intp_data_new.resize(XDoF());
		}

		//number of cols
		virtual int XDoF() const {
			return fine_cell_type->grid.Counts().prod();
		}

		//number of rows
		virtual int YDoF() const {
			return coarse_cell_type->grid.Counts().prod();
		}

		//input p, get Ap
		virtual void Apply(ArrayDv<T>& coarse_data, const ArrayDv<T>& fine_data) {
			Base::Memory_Check(coarse_data, fine_data, "Restrictor::Apply error: not enough space");
			
			T* intp_ptr_old = ArrayFunc::Data(intp_data_old);
			T* intp_ptr_new = ArrayFunc::Data(intp_data_new);
			const T* fine_data_ptr = ArrayFunc::Data(fine_data);
			const unsigned char* fine_cell_type_ptr = fine_cell_type->Data_Ptr();
			
			cudaMemcpy(intp_ptr_old, fine_data_ptr, sizeof(T) * fine_data.size(), cudaMemcpyDeviceToDevice);

			auto set_fixed = [=]__device__(T & a, const  unsigned char& type) { if (type == 1 || type == 2)a = 0; };
			GPUFunc::Cwise_Mapping_Wrapper(intp_ptr_old, fine_cell_type_ptr, set_fixed, fine_data.size());

			for (int axis = 0; axis < d; axis++) {
				fine_cell_type->Exec_Kernel(&Restrictor_Intp_Axis_Kernel<T, d>, axis, fine_cell_type->grid, intp_ptr_new, intp_ptr_old);
				std::swap(intp_ptr_old, intp_ptr_new);
			}
			coarse_cell_type->Exec_Kernel(&Restrictor_Intp_Coarser_Kernel<T, d>, coarse_cell_type->grid, ArrayFunc::Data(coarse_data), fine_cell_type->grid, intp_ptr_old);
		
			T* coarse_data_ptr = ArrayFunc::Data(coarse_data);
			const unsigned char* coarse_cell_type_ptr = coarse_cell_type->Data_Ptr();
			GPUFunc::Cwise_Mapping_Wrapper(coarse_data_ptr, coarse_cell_type_ptr, set_fixed, coarse_data.size());
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
			VectorDi fine_coord = coarse_coord * 2 + MathFunc::Vi<d>(dx, dy, dz);
			if (fine_grid.Valid(fine_coord)) sum += fine_data[fine_grid.Index(fine_coord)];
		}
		coarse_data[coarse_grid.Index(coarse_coord)] = sum / bc_num * 4;
	}

	template<class T, int d>
	class RestrictorSum : public LinearMapping<T> {
	public:
		Grid<d> coarse_grid, fine_grid;

		RestrictorSum() {}
		RestrictorSum(const FieldDv<unsigned char, d>& coarse_cell_type, const FieldDv<unsigned char, d>& fine_cell_type) {
			Init(coarse_cell_type, fine_cell_type);
		}
		void Init(const FieldDv<unsigned char, d>& coarse_cell_type, const FieldDv<unsigned char, d>& fine_cell_type) {
			coarse_grid = coarse_cell_type.grid, fine_grid = fine_cell_type.grid;
		}

		//number of cols
		virtual int XDoF() const { return fine_grid.Valid_Size(); }
		//number of rows
		virtual int YDoF() const { return coarse_grid.Valid_Size(); }

		virtual void Apply(ArrayDv<T>& coarse_data, const ArrayDv<T>& fine_data) {
			Memory_Check(coarse_data, fine_data, "RestrictorSum::Apply error: not enough space");
			T* coarse_ptr = ArrayFunc::Data(coarse_data);
			const T* fine_ptr = ArrayFunc::Data(fine_data);
			coarse_grid.Exec_Kernel(&Restrictor_Sum_Coarse_Kernel<T, d>, coarse_grid, coarse_ptr, fine_grid, fine_ptr);
			checkCudaErrors(cudaGetLastError());
		}
	};
}