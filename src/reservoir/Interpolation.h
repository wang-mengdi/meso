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
	
	template<class T, int d, DataHolder side> class FaceField;

	class PointIntpLinear {
	public:
		template<class T, int d>
		static T __host__ __device__ Value(const Grid<d> grid, const T* data, const Vector<int, d> coord, const Vector<real, d> frac) {
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
			else Assert("PointIntpLinear::Value error: dimension must be 2 or 3");
		}
	};

	class PointIntpLinearPadding0 {
	public:
		template<class T, int d>
		static T __host__ __device__ Value(const Grid<d> grid, const T* data, const Vector<int, d> coord, const Vector<real, d> frac) {
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
	};

	class PointIntpLinearClamp {
	public:
		template<class T, int d>
		static T __host__ __device__ Value(const Grid<d> grid, const T* data, const Vector<int, d> coord, const Vector<real, d> frac) {
			static constexpr T padding_val = 0;
			//considering invalid datas as 0
			Typedef_VectorD(d);
			static constexpr int dx[8] = { 0,1,0,1,0,1,0,1 };
			static constexpr int dy[8] = { 0,0,1,1,0,0,1,1 };
			static constexpr int dz[8] = { 0,0,0,0,1,1,1,1 };
			//clamping
			for (int axis = 0; axis < d; axis++) {
				if (coord[axis] < 0) coord[axis] = 0, frac[axis] = 0;
				if (coord[axis] > grid.counts[axis] - 2) coord[axis] = grid.counts[axis] - 2, frac[axis] = 1;
			}
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
			else Assert("PointIntpLinearClamp::Value error: dimension must be 2 or 3");
		}
	};

	template<class PointIntp>
	class Interpolation {
	public:
		Interpolation() {}
		template<class T, int d>
		static T __host__ __device__ Value(const Grid<d> grid, const T* data, const Vector<int, d> coord, const Vector<real, d> frac) {
			return PointIntp::Value(grid, data, coord, frac);
		}
		template<class T, int d>
		static T __host__ __device__ Value(const Grid<d> grid, const T* data, const Vector<real, d> pos) {
			Vector<int, d> node;
			Vector<real, d> frac;
			grid.Get_Fraction(pos, node, frac);
			return PointIntp::Value(grid, data, node, frac);
		}
		template<class T, int d, DataHolder side>
		static T __host__ __device__ Value(const Field<T, d, side>& F, const Vector<real, d> pos) {
			const T* data_ptr = F.Data_Ptr();
			Vector<int, d> node; Vector<real, d> frac;
			F.grid.Get_Fraction(pos, node, frac);
			return PointIntp::Value(F.grid, data_ptr, node, frac);
		}
		template<class T, int d>
		static Vector<T, d> __host__ __device__ Face_Vector(const Grid<d> g0, const T* v0, const Grid<d> g1, const T* v1, const Grid<d> g2, const T* v2, const Vector<real, d> pos) {
			Typedef_VectorD(d);
			Vector<T, d> ret;
			VectorDi node; VectorD frac;
			//x
			{
				g0.Get_Fraction(pos, node, frac);
				ret[0] = PointIntp::Value(g0, v0, node, frac);
			}
			//y
			if constexpr (d >= 2) {
				g1.Get_Fraction(pos, node, frac);
				ret[1] = PointIntp::Value(g1, v1, node, frac);
			}
			//z
			if constexpr (d >= 3) {
				g2.Get_Fraction(pos, node, frac);
				ret[2] = PointIntp::Value(g2, v2, node, frac);
			}
			return ret;
		}
		template<class T, int d, DataHolder side>
		static Vector<T, d> __host__ __device__ Face_Vector(const FaceField<T, d, side>& vector_field, const Vector<real, d> pos) {
			const auto& grid = vector_field.grid;
			Grid<d> g0 = grid.Face_Grid(0), g1 = grid.Face_Grid(1), g2 = grid.Face_Grid(2);
			const T* v0 = vector_field.Data_Ptr(0), * v1 = vector_field.Data_Ptr(1), * v2 = vector_field.Data_Ptr(2);
			return Interpolation<PointIntp>::Face_Vector(g0, v0, g1, v1, g2, v2, pos);
		}
	}; 

	using IntpLinear = Interpolation<PointIntpLinear>;
	using IntpLinearPadding0 = Interpolation<PointIntpLinearPadding0>;
	using IntpLinearClamp = Interpolation<PointIntpLinearClamp>;
}