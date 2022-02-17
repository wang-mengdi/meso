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
		FieldDv<T, d> vol;
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

		virtual int cols() const { return dof; }//number of cols

		virtual int Y_DoF() const { return dof; }//number of rows

		//input p, get Ap
		virtual void Apply(ArrayDv<T>& Ap, const ArrayDv<T>& p) {
			//all data is saved in d_temp. as face data and cell data(p)
			//int face_x_off = 0;
			//int face_y_off = face_x_off + grid.face_size(0);
			//int cell_off = face_y_off + grid.face_size(1);

			//copy p to temp
			ArrayFunc::Copy(temp_cell.data, p);

			auto fix = [=] __device__ (T v, bool fixed) ->T {return fixed ? 0 : v; };

			thrust::transform(temp_cell.data.begin(), temp_cell.data.end(), fixed.data.begin(), temp_cell.data.begin(), fix);
			//BinaryOp<T, bool, T> my_fix(fix);
			//ArrayFunc::Binary_Transform(fix, temp_cell.data, fixed.data, temp_cell.data);



			//cudaMemset(Ap, 0, sizeof(Scalar) * dof);

			//auto fix = [=] __device__(Scalar & tv, bool tfixed) { if (tfixed) tv = (Scalar)0; };
			//auto multi = [=] __device__(Scalar & tv, Scalar tvol) { tv *= tvol; };
			//auto neg = [=]__device__(Scalar & tv) { tv = -tv; };



			////save p in device
			//cudaMemcpy(d_temp + cell_off, p, sizeof(Scalar) * grid.cell_size(), cudaMemcpyDeviceToDevice);
			//cwise_mapping_wrapper(d_temp + cell_off, descr->d_fixed, fix, grid.cell_size());

			//grid2DOperator::Cod0Mapping(grid, d_temp + face_x_off, d_temp + face_y_off, d_temp + cell_off);

			//cwise_mapping_wrapper(d_temp + face_y_off, neg, grid.face_size(1));
			//cwise_mapping_wrapper(d_temp + face_x_off, descr->d_vol + face_x_off, multi, grid.face_size(0));
			//cwise_mapping_wrapper(d_temp + face_y_off, descr->d_vol + face_y_off, multi, grid.face_size(1));

			//grid2DOperator::D1Mapping(grid, Ap, d_temp + face_x_off, d_temp + face_y_off);

			//cwise_mapping_wrapper(Ap, descr->d_fixed, fix, grid.cell_size());

			//cwise_mapping_wrapper(Ap, neg, grid.cell_size());

			//auto cond_add = [=]__device__(Scalar & tv1, Scalar tv2, bool tfixed) { if (tfixed) tv1 += tv2; };
			//cwise_mapping_wrapper(Ap, p, descr->d_fixed, cond_add, grid.cell_size());
		}
	};

}