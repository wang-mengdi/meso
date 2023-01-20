//////////////////////////////////////////////////////////////////////////
// Poisson Linear Mapping with celltype and vol
// Copyright (c) (2022-), Zangyueyang Xian, Mengdi Wang, Yuchen Sun
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "LinearMapping.h"
#include "Field.h"
#include "FaceField.h"
#include "ExteriorDerivative.h"
#include "AuxFunc.h"
#include "BoundaryCondition.h"
using namespace thrust::placeholders;

namespace Meso {
	//Negative Poisson mapping -lap(p), except some masked points
	//Masked cells will be viewed as 0 in poisson mapping
	//Which means adjacent faces of masked cells will have volume 0
	template<class T, int d, class ExteriorDerivative=ExteriorDerivativePadding0>
	class MaskedPoissonMapping : public LinearMapping<T> {
		using Base=LinearMapping<T>;
		//Ap=-lap(p)
		Typedef_VectorD(d);
	public:
		int dof;
		FieldDv<unsigned char, d> cell_type;
		FaceFieldDv<T, d> vol;
		ArrayDv<bool> is_boundary_tile;
		ArrayDv<int> boundary_tiles;

		FaceFieldDv<T, d> temp_face;
		
		MaskedPoissonMapping() {}
		MaskedPoissonMapping(const Grid<d> grid) { Allocate_Memory(grid); }

		void Allocate_Memory(const Grid<d> grid) {
			Assert(grid.Is_Unpadded(), "MaskedPoissonMapping: invalid grid {}, padding not allowed", grid);
			dof = grid.Memory_Size();
			cell_type.Init(grid);
			vol.Init(grid);
			temp_face.Init(grid);
		}

		void Init(const Grid<d> grid) {
			Allocate_Memory(grid);
		}

		void Init(const Field<unsigned char, d>& _cell_type, const FaceField<T, d>& _vol) {
			Assert(_cell_type.grid.Indexer() == _vol.grid.Indexer(), "MaskedPoissonMapping::Init error: _cell_type grid {} unequal to _vol grid {}", _cell_type.grid, _vol.grid);
			Allocate_Memory(_cell_type.grid);
			cell_type.Deep_Copy(_cell_type);
			vol.Deep_Copy(_vol);
			
		}

		const Grid<d>& Grid(void) const {
			return vol.grid;
		}

		virtual int XDoF() const { return dof; }//number of cols

		virtual int YDoF() const { return dof; }//number of rows

		//input p, get Ap
		virtual void Apply(ArrayDv<T>& Ap, const ArrayDv<T>& p) {
			Base::Memory_Check(Ap, p, "PoissonMapping::Apply error: not enough space");

			auto fix = [=] __device__(T & tv,  unsigned char type) { if (type == 1 || type == 2) tv = (T)0; };
			auto multi = [=] __device__(T & tv, T tvol) { tv *= tvol; };
			auto neg = [=]__device__(T & tv) { tv = -tv; };
			auto cond_set = [=]__device__(T & tv1, T tv2, unsigned char type) { if (type == 1 || type == 2) tv1 = tv2; };

			Meso::Grid<d> grid = Grid();
			T* Ap_ptr = thrust::raw_pointer_cast(Ap.data());
			const T* p_ptr = thrust::raw_pointer_cast(p.data());

			//temp_cell=p, set to 0 for dirichlet and neumann cell
			cudaMemcpy(Ap_ptr, p_ptr, sizeof(T)* dof, cudaMemcpyDeviceToDevice);
			GPUFunc::Cwise_Mapping_Wrapper(Ap_ptr, cell_type.Data_Ptr(), fix, dof);

			//temp_face = grad(temp_cell) *. vol
			//d(p) ----- 1-form
			for (int axis = 0; axis < d; axis++)
				cudaMemset(temp_face.Data_Ptr(axis), 0, sizeof(T) * grid.Face_Grid(axis).Memory_Size());
			if constexpr (d == 2) grid.Exec_Kernel(&D_CoCell_Kernel2_Padding0<T>, grid, temp_face.Data_Ptr(0), temp_face.Data_Ptr(1), Ap_ptr);
			else if constexpr (d == 3) grid.Exec_Kernel(&D_CoCell_Kernel3_Padding0<T>, grid, temp_face.Data_Ptr(0), temp_face.Data_Ptr(1), temp_face.Data_Ptr(2), Ap_ptr);
			
			//d(p) *. vol ----- 1-form
			for (int axis = 0; axis < d; axis++)
				GPUFunc::Cwise_Mapping_Wrapper(temp_face.Data_Ptr(axis), vol.Data_Ptr(axis), multi, grid.Face_Grid(axis).Memory_Size());
			
			//div(temp_face)
			//d*d(p) *. vol ----- 3-form
			if constexpr (d == 2) grid.Exec_Kernel(D_Face_Kernel2_Padding0<T>, grid, Ap_ptr, temp_face.Data_Ptr(0), temp_face.Data_Ptr(1));
			else if constexpr (d == 3)grid.Exec_Kernel(D_Face_Kernel3_Padding0<T>, grid, Ap_ptr, temp_face.Data_Ptr(0), temp_face.Data_Ptr(1), temp_face.Data_Ptr(2));
			
			//neg div
			GPUFunc::Cwise_Mapping_Wrapper(Ap_ptr, neg, dof);
			
			//transfer dirichlet and neumann cell data to Ap
			GPUFunc::Cwise_Mapping_Wrapper(Ap_ptr, p_ptr, cell_type.Data_Ptr(), cond_set, dof);
		}


	};

}
