//////////////////////////////////////////////////////////////////////////
// Poisson Linear Mapping with Fixed Mask
// Copyright (c) (2022-), Zangyueyang Xian, Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "LinearMapping.h"
#include "Field.h"
#include "FaceField.h"
#include "LambdaHelper.h"
#include "DifferentialExteriorCalculus.h"
#include "AuxFunc.h"
#include "BoundaryCondition.h"
using namespace thrust::placeholders;

namespace Meso {
	template<class T, int d>
	class MaskedPoissonMapping : public LinearMapping<T> {
		//Ap=-lap(p)
		Typedef_VectorD(d);
	public:
		int dof;
		FaceFieldDv<T, d> vol;
		//BoundaryConditionDirect<FieldDv<bool, d>> cell_bc;
		FieldDv<bool, d> fixed;

		FieldDv<T, d> temp_cell;
		FaceFieldDv<T, d> temp_face;
		
		MaskedPoissonMapping() {}
		MaskedPoissonMapping(const Grid<d>& grid) { Allocate_Memory(grid); }

		void Allocate_Memory(const Grid<d>& grid) {
			dof = grid.DoF();
			vol.Init(grid);
			fixed.Init(grid);
			temp_cell.Init(grid);
			temp_face.Init(grid);
		}
		void Init(const Grid<d>& grid) {
			Allocate_Memory(grid);
		}
		void Init(const Grid<d>& grid, const FaceField<T, d>& _vol, const Field<bool, d>& _fixed) {
			Allocate_Memory(grid);
			vol.Deep_Copy(_vol);
			//cell_bc.Init(_fixed);
			fixed.Deep_Copy(_fixed);
		}
		//template<class IFFunc, class CFunc>
		//void Init(const Grid<d>& grid, IFFunc vol_func, CFunc is_unknown_func) {
		//	Allocate_Memory(grid);
		//	vol.Calc_Faces(vol_func);
		//	fixed.Calc_Cells(
		//		[=](const VectorDi& cell)->bool {return !is_unknown_func(cell); }
		//	);
		//}

		const Grid<d>& Grid(void) const {
			return vol.grid;
		}

		virtual int XDof() const { return dof; }//number of cols

		virtual int YDof() const { return dof; }//number of rows

		//input p, get Ap
		virtual void Apply(ArrayDv<T>& Ap, const ArrayDv<T>& p) {
			//Poisson mapping is *d*d in DEC(Discrete Exterior Calculus) theory
			//so lap(p)=*d*d(p)

			Memory_Check(Ap, p, "PoissonMapping::Apply error: not enough space");

			//temp_cell=p, set to 0 if fixed			
			//ArrayFunc::Copy(temp_cell.Data(), p);
			temp_cell.Fill(0);
			ArrayFunc::Copy_UnMasked(temp_cell.Data(), p, fixed.Data());


			//temp_face = grad(temp_cell) *. vol
			//d(p) ----- 1-form
			Exterior_Derivative(temp_face, temp_cell);
			//d(p) *. vol ----- 1-form
			temp_face *= vol;

			//Hodge star is identity here
			//*d(p) *. vol ----- 2-form

			//temp_cell = div(temp_face)
			//d*d(p) *. vol ----- 3-form
			Exterior_Derivative(temp_cell, temp_face);
			//temp_cell *= -1;

			//Hodge star is identity here
			//*d*d(p) *. vol ----- 0-form

			//transfer data to Ap
			ArrayFunc::Copy(Ap, temp_cell.Data());
			ArrayFunc::Copy_Masked(Ap, p, fixed.Data());

			checkCudaErrors(cudaGetLastError());
		}
	};

}