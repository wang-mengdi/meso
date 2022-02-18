//////////////////////////////////////////////////////////////////////////
// Poisson Linear Mapping
// Copyright (c) (2018-), Zangyueyang Xian, Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "LinearMapping.h"
#include "Field.h"
#include "FaceField.h"
#include "LambdaHelper.h"
#include "DifferentialGeometry.h"
using namespace thrust::placeholders;

namespace Meso {

	template<class T1, class T2, class T3>
	class BinaryOp {
		std::function<T3(T1, T2)> f;
		template<class TTFuncT>
		BinaryOp(TTFuncT _f) {
			f = _f;
		}
		__host__ __device__ T3 operator () (const T1 a, const T2 b) {
			return f(t1, t2);
		}
	};

	template<class T, int d>
	class PoissonMapping : public LinearMapping<T> {
		Typedef_VectorD(d);
	public:
		int dof;
		FaceFieldDv<T, d> vol;
		FieldDv<bool, d> fixed;

		FieldDv<T, d> temp_cell;
		FaceFieldDv<T, d> temp_face;

		void Init(const Grid<d, GridType::CELL>& grid, IFFunc<T, d> vol_func, CFunc<T, d> is_unknown_func) {
			dof = grid.DoF();
			vol.Calc_Each(vol_func);
			fixed.Calc_Each(
				[=](const VectorDi& cell)->bool {return !is_unknown_func(cell); }
			);
			temp_cell.Init(grid);
			temp_face.Init(grid);
		}

		virtual int XDof() const { return dof; }//number of cols

		virtual int YDof() const { return dof; }//number of rows

		//input p, get Ap
		virtual void Apply(ArrayDv<T>& Ap, const ArrayDv<T>& p) {
			//temp_cell=p, set to 0 if fixed
			auto identity_except_fixed = [=] __device__(T v, bool fixed) ->T { return fixed ? 0 : v; };
			ArrayFunc::Binary_Transform(p, fixed.data, identity_except_fixed, temp_cell.data);

			//temp_face = grad(temp_cell) *. vol
			D_CoCell_Mapping(temp_cell, temp_face);
			ArrayFunc::Unary_Transform(temp_face.face_data[1], thrust::negate<T>(), temp_face.face_data[1]);
			ArrayFunc::Binary_Transform(temp_face.face_data[0], vol.face_data[0], thrust::multiplies<T>(), temp_face.face_data[0]);
			ArrayFunc::Binary_Transform(temp_face.face_data[1], vol.face_data[1], thrust::multiplies<T>(), temp_face.face_data[1]);

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