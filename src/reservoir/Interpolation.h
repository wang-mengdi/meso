//////////////////////////////////////////////////////////////////////////
// Interpolation on grid
// Copyright (c) (2022-), Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once

#include "AuxFunc.h"
#include "Grid.h"
#include "Field.h"
#include "FaceField.h"

namespace Meso {
	//template<class T, int d, GridType gtype=CENTER>
	//class Interpolation {
	//public:
	//	virtual static T __host__ __device__ Point_Interpolate(const Grid<d, gtype> grid, const T* data, const Vector<int, d> coord, const Vector<real, d> frac) = 0;

	//};

	namespace Interpolation {
		template<class T, int d, GridType gtype>
		T __host__ __device__ Linear_Intp(const Grid<d, gtype> grid, const T* data, const Vector<int, d> coord, const Vector<real, d> frac) {
			Typedef_VectorD(d);
			static constexpr int dx[8] = { 0,1,0,1,0,1,0,1 };
			static constexpr int dy[8] = { 0,0,1,1,0,0,1,1 };
			static constexpr int dz[8] = { 0,0,0,0,1,1,1,1 };
			if constexpr (d == 2) {
				real w[2][2] = { {1.0 - frac[0],frac[0]},{1.0 - frac[1],frac[1]} };
				T intp_value = 0;
				for (int s = 0; s < 4; s++) {
					int d0 = dx[s], d1 = dy[s];
					int idx = grid.Index(coord[0] + d0, coord[1] + d1);
					intp_value += w[0][d0] * w[1][d1] * data[idx];
				}
				return intp_value;
			}
			else if constexpr (d == 3) {
				real w[3][2] = { {1.0 - frac[0],frac[0]},{1.0 - frac[1],frac[1]} ,{1.0 - frac[2],frac[2]} };
				T intp_value = 0;
				for (int s = 0; s < 8; s++) {
					int d0 = dx[s], d1 = dy[s], d2 = dz[s];
					int idx = grid.Index(coord[0] + d0, coord[1] + d1, coord[2] + d2);
					intp_value += w[0][d0] * w[1][d1] * w[2][d2] * data[idx];
				}
				return intp_value;
			}
			else Assert("Interpolation:Linear_Intp error: dimension must be 2 or 3");
		}
		template<class T,int d, DataHolder side>
		T __host__ __device__ Linear_Intp(const Field<T, d, side>& F, const Vector<real, d> pos) {
			const T* data_ptr = F.Data_Ptr();
			Vector<int, d> node; Vector<real, d> frac;
			F.grid.Get_Fraction(pos, node, frac);
			return Linear_Intp(F.grid, data_ptr, node, frac);
		}

		template<class T, int d, GridType gtype>
		T __host__ __device__ Linear_Intp_Padding0(const Grid<d, gtype> grid, const T* data, const Vector<int, d> coord, const Vector<real, d> frac) {
			static constexpr T padding_val = 0;
			//considering invalid datas as 0
			Typedef_VectorD(d);
			static constexpr int dx[8] = { 0,1,0,1,0,1,0,1 };
			static constexpr int dy[8] = { 0,0,1,1,0,0,1,1 };
			static constexpr int dz[8] = { 0,0,0,0,1,1,1,1 };
			if constexpr (d == 2) {
				real w[2][2] = { {1.0 - frac[0],frac[0]},{1.0 - frac[1],frac[1]} };
				T intp_value = 0;
				for (int s = 0; s < 4; s++) {
					int d0 = dx[s], d1 = dy[s];
					T val = grid.Valid(coord[0] + d0, coord[1] + d1) ? data[grid.Index(coord[0] + d0, coord[1] + d1)] : padding_val;
					intp_value += w[0][d0] * w[1][d1] * val;
				}
				return intp_value;
			}
			else if constexpr (d == 3) {
				real w[3][2] = { {1.0 - frac[0],frac[0]},{1.0 - frac[1],frac[1]} ,{1.0 - frac[2],frac[2]} };
				T intp_value = 0;
				for (int s = 0; s < 8; s++) {
					int d0 = dx[s], d1 = dy[s], d2 = dz[s];
					T val = grid.Valid(coord[0] + d0, coord[1] + d1, coord[2] + d2) ? data[grid.Index(coord[0] + d0, coord[1] + d1, coord[2] + d2)] : padding_val;
					intp_value += w[0][d0] * w[1][d1] * w[2][d2] * val;
				}
				return intp_value;
			}
			else Assert("Interpolation:Linear_Intp error: dimension must be 2 or 3");
		}
		template<class T, int d, GridType gtype, DataHolder side>
		Vector<T, d> Linear_Intp_Vector_Padding0(const FaceField<T, d, side>& vector_field, const Vector<real, d> pos) {
			Typedef_VectorD(d);
			Vector<T, d> ret;
			for (int axis = 0; axis < d; axis++) {
				VectorDi node; VectorD frac;
				vector_field.grid.Get_Fraction(pos, node, frac);
				ret[axis] = Linear_Intp_Padding0(vector_field.grid, vector_field.Data_Ptr(axis), node, frac);
			}
			return ret;
		}
	}
}