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
		
		void Allocate_Memory(const Grid<d, CELL>& grid) {
			dof = grid.DoF();
			vol.Init(grid);
			fixed.Init(grid);
			temp_cell.Init(grid);
			temp_face.Init(grid);
		}
		void Init(const Grid<d, CELL>& grid, const FaceField<T, d>& _vol, const Field<bool, d>& _fixed) {
			Allocate_Memory(grid);
			vol.Copy(_vol);
			fixed.Copy(_fixed);
		}
		template<class IFFunc, class CFunc>
		void Init(const Grid<d, GridType::CELL>& grid, IFFunc vol_func, CFunc is_unknown_func) {
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
			//temp_cell=p, set to 0 if fixed
			auto identity_except_fixed = [=] __device__(T v, bool fixed) ->T { return fixed ? 0 : v; };
			ArrayFunc::Binary_Transform(p, fixed.data, identity_except_fixed, temp_cell.data);

			Info("temp_cell [7,0] {}", temp_cell.data[temp_cell.grid.Index(Vector2i(7, 0))]);
			Info("temp_cell [8,0] {}", temp_cell.data[temp_cell.grid.Index(Vector2i(8, 0))]);

			Info("before cocell_mapping 0 [8 0] {}", temp_face.face_data[0][vol.grid.Face_Index(0, Vector2i(8, 0))]);

			//temp_face = grad(temp_cell) *. vol
			D_CoCell_Mapping(temp_cell, temp_face);

			Info("after cocell_mapping 0 [8 0] {}", temp_face.face_data[0][vol.grid.Face_Index(0, Vector2i(8, 0))]);

			ArrayFunc::Unary_Transform(temp_face.face_data[1], thrust::negate<T>(), temp_face.face_data[1]);
			ArrayFunc::Binary_Transform(temp_face.face_data[0], vol.face_data[0], thrust::multiplies<T>(), temp_face.face_data[0]);
			ArrayFunc::Binary_Transform(temp_face.face_data[1], vol.face_data[1], thrust::multiplies<T>(), temp_face.face_data[1]);

			if constexpr (d == 2) {
				Info("7 neighbors temp_face:");
				//Info("0 [7 0] {}", temp_face(0, Vector2i(7, 0)));
				Info("0 [7 0] {}", temp_face.face_data[0][vol.grid.Face_Index(0, Vector2i(7, 0))]);
				Info("0 [8 0] {}", temp_face.face_data[0][vol.grid.Face_Index(0, Vector2i(8, 0))]);
				Info("1 [7 0] {}", temp_face.face_data[1][vol.grid.Face_Index(1, Vector2i(7, 0))]);
				Info("1 [7 1] {}", temp_face.face_data[1][vol.grid.Face_Index(1, Vector2i(7, 1))]);
				Info("vol 0 [8 0] {}", vol.face_data[0][vol.grid.Face_Index(0, Vector2i(8, 0))]);
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