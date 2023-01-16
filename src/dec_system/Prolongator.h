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
	__global__ void Prolongator_Intp_Kernel(const GridIndexer<d> fine_grid, T* fine_data, const GridIndexer<d> coarse_grid, const T* coarse_data) {
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
		fine_data[fine_grid.Index(fine_coord)] = IntpLinearPadding0::Value(coarse_grid, coarse_data, coarse_coord, coarse_frac);
		//fine_data[fine_grid.Index(fine_coord)] = Interpolation::Linear_Intp_Padding0(coarse_grid, coarse_data, coarse_coord, coarse_frac);
	}

	template<class T, int d>
	class ProlongatorIntp : public LinearMapping<T> {
		Typedef_VectorD(d);
	public:
		using Base=LinearMapping<T>;
		//GridIndexer<d> fine_grid, coarse_grid;
		const FieldDv<CellType, d>* fine_cell_type;
		const FieldDv<CellType, d>* coarse_cell_type;
		ArrayDv<T> temp_data;

		ProlongatorIntp() {}
		ProlongatorIntp(const FieldDv<CellType, d>& _fine_cell_type, const FieldDv<CellType, d>& _coarse_cell_type) { Init(_fine_cell_type, _coarse_cell_type); }

		void Init(const FieldDv<CellType, d>& _fine_cell_type, const FieldDv<CellType, d>& _coarse_cell_type) {
			Assert(_fine_cell_type.grid.Is_Unpadded(), "ProlongatorIntp: _fine {} invalid, must be unpadded", _fine_cell_type.grid);
			Assert(_coarse_cell_type.grid.Is_Unpadded(), "ProlongatorIntp: _coarse {} invalid, must be unpadded", _coarse_cell_type.grid);
			fine_cell_type = &_fine_cell_type;
			coarse_cell_type = &_coarse_cell_type;
			temp_data.resize(XDoF());
		}
		//number of cols
		virtual int XDoF() const {
			return coarse_cell_type->grid.Counts().prod();
		}

		//number of rows
		virtual int YDoF() const {
			return fine_cell_type->grid.Counts().prod();
		}

		//input p, get Ap
		virtual void Apply(ArrayDv<T>& fine, const ArrayDv<T>& coarse) {
			Base::Memory_Check(fine, coarse, "Prolongator::Apply error: not enough memory");

			T* temp_data_ptr = ArrayFunc::Data(temp_data);
			T* fine_ptr = ArrayFunc::Data(fine);
			const T* coarse_ptr = ArrayFunc::Data(coarse);
			const CellType* coarse_cell_type_ptr = coarse_cell_type->Data_Ptr();
			const CellType* fine_cell_type_ptr = fine_cell_type->Data_Ptr();

			cudaMemcpy(temp_data_ptr, coarse_ptr, sizeof(T) * coarse.size(), cudaMemcpyDeviceToDevice);
			
			auto set_fixed = [=]__device__(T & a, const  CellType& type) { if (type == AIR || type == SOLID) a = 0; };
			GPUFunc::Cwise_Mapping_Wrapper(temp_data_ptr, coarse_cell_type_ptr, set_fixed, temp_data.size());

			fine_cell_type->Exec_Kernel(&Prolongator_Intp_Kernel<T, d>, fine_cell_type->grid, fine_ptr, coarse_cell_type->grid, coarse_ptr);
			
			GPUFunc::Cwise_Mapping_Wrapper(fine_ptr, fine_cell_type_ptr, set_fixed, fine.size());
			
			checkCudaErrors(cudaGetLastError());
		}
	};

	template<class T, int d>
	__global__ void Prolongator_Sum_Kernel(const Grid<d> fine_grid, T* fine_data, const Grid<d> coarse_grid, const T* coarse_data) {
		Typedef_VectorD(d);
		VectorDi fine_coord = GPUFunc::Thread_Coord<d>(blockIdx, threadIdx);
		VectorDi coarse_coord = fine_coord / 2;
		real val = 0;
		if (coarse_grid.Valid(coarse_coord)) val = coarse_data[coarse_grid.Index(coarse_coord)];
		fine_data[fine_grid.Index(fine_coord)] = val;
	}

	template<class T, int d>
	class ProlongatorSum : public LinearMapping<T> {
	public:
		Grid<d> fine_grid, coarse_grid;

		ProlongatorSum() {}
		ProlongatorSum(const Grid<d> _fine, const Grid<d> _coarse) { Init(_fine, _coarse); }
		void Init(const Grid<d> _fine, const Grid<d> _coarse) {
			fine_grid = _fine;
			coarse_grid = _coarse;
		}

		//number of cols
		virtual int XDoF() const { return coarse_grid.DoF(); }
		//number of rows
		virtual int YDoF() const { return fine_grid.DoF(); }
		//input coarse_data, output fine_data
		virtual void Apply(ArrayDv<T>& fine_data, const ArrayDv<T>& coarse_data) {
			Memory_Check(fine_data, coarse_data, "ProlongatorSum::Apply error: not enough memory");
			T* fine_ptr = ArrayFunc::Data(fine_data);
			const T* coarse_ptr = ArrayFunc::Data(coarse_data);
			fine_grid.Exec_Kernel(&Prolongator_Sum_Kernel<T, d>, fine_grid, fine_ptr, coarse_grid, coarse_ptr);
			checkCudaErrors(cudaGetLastError());
		}
	};
}
