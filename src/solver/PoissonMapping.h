//////////////////////////////////////////////////////////////////////////
// Poisson Linear Mapping
// Copyright (c) (2022-), Zangyueyang Xian, Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "LinearMapping.h"
#include "Field.h"
#include "FaceField.h"
#include "LambdaHelper.h"
#include "DifferentialGeometry.h"
#include "AuxFunc.h"
using namespace thrust::placeholders;

namespace Meso {
	template<class T, int d>
	class PoissonMapping : public LinearMapping<T> {
		//Ap=-lap(p)
		Typedef_VectorD(d);
	public:
		int dof;
		FaceFieldDv<T, d> vol;
		FieldDv<bool, d> fixed;

		FieldDv<T, d> temp_cell;
		FaceFieldDv<T, d> temp_face;
		
		void Allocate_Memory(const Grid<d, CENTER>& grid) {
			dof = grid.DoF();
			vol.Init(grid);
			fixed.Init(grid);
			temp_cell.Init(grid);
			temp_face.Init(grid);
		}
		void Init(const Grid<d, CENTER>& grid, const FaceField<T, d>& _vol, const Field<bool, d>& _fixed) {
			Allocate_Memory(grid);
			vol.Copy(_vol);
			fixed.Copy(_fixed);
		}
		template<class IFFunc, class CFunc>
		void Init(const Grid<d, CENTER>& grid, IFFunc vol_func, CFunc is_unknown_func) {
			Allocate_Memory(grid);
			vol.Calc_Faces(vol_func);
			fixed.Calc_Cells(
				[=](const VectorDi& cell)->bool {return !is_unknown_func(cell); }
			);
		}

		virtual int XDof() const { return dof; }//number of cols

		virtual int YDof() const { return dof; }//number of rows

		//input p, get Ap
		virtual void Apply(ArrayDv<T>& Ap, const ArrayDv<T>& p) {
			Assert(p.size() == dof, "PoissonMapping: p.size() not equal to dof");
			Assert(Ap.size() == dof, "PoissonMapping: Ap.size() not equal to dof");

			//temp_cell=p, set to 0 if fixed
			auto identity_except_fixed = [=] __device__(T v, bool fixed) ->T { return fixed ? 0 : v; };
			ArrayFunc::Binary_Transform(p, fixed.data, identity_except_fixed, temp_cell.data);

			//temp_face = grad(temp_cell) *. vol
			D_CoCell_Mapping(temp_cell, temp_face);
			for (int axis = 0; axis < d; axis++) {
				ArrayFunc::Multiply(temp_face.face_data[axis], vol.face_data[axis]);
			}

			//temp_cell = -div(temp_face)
			D_Face_Mapping(temp_face, temp_cell);
			ArrayFunc::Unary_Transform(temp_cell.data, thrust::negate<T>(), temp_cell.data);
			//Ap=temp_cell, set to 0 if fixed
			ArrayFunc::Binary_Transform(temp_cell.data, fixed.data, identity_except_fixed, Ap);

			//if fixed, add p back
			thrust::transform_if(
				Ap.begin(),//first1
				Ap.end(),//last1
				p.begin(),//first2
				fixed.data.begin(),//stencil
				Ap.begin(),//result
				_1 + _2,//binary op
				thrust::identity<bool>()//pred
			);
		}
	};

}