//////////////////////////////////////////////////////////////////////////
// Interpolation on grid, but with extra data reference
// Copyright (c) (2022-), Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once

#include "AuxFunc.h"
#include "Field.h"
#include "FaceField.h"

namespace Meso {
	template<class T, int d, DataHolder side>
	class Interpolator {
		Typedef_VectorD(d);
	public:
		__host__ __device__ virtual T Value(const GridIndexer<d> grid, const T* data, const Vector<int, d> coord, const Vector<real, d> frac) = 0;
		__host__ __device__ T Value1(Field<T, d, side>&& F, const Vector<real, d> pos) {
			const T* data_ptr = F.Data_Ptr();
			Vector<int, d> node; Vector<real, d> frac;
			F.grid.Get_Fraction(pos, node, frac);
			return Value(F.grid, data_ptr, node, frac);
		}
	};

	template<class T, int d, DataHolder side>
	class InterpolatorLinearWithMask : public Interpolator<T, d, side> {
		Typedef_VectorD(d);
	public:
		//only interpolate on cells which are valid and valid_mask==true 
		Field<bool, d, side> valid_mask;

		void Init_Shallow(const Field<bool, d, side>& _valid_mask) {
			valid_mask.Shallow_Copy(_valid_mask);
		}

		__host__ __device__ T Value(const GridIndexer<d> grid, const T* data, const Vector<int, d> coord, const Vector<real, d> frac) {
			static constexpr T padding_val = 0;
			//considering invalid datas as 0
			static constexpr int dx[8] = { 0,1,0,1,0,1,0,1 };
			static constexpr int dy[8] = { 0,0,1,1,0,0,1,1 };
			static constexpr int dz[8] = { 0,0,0,0,1,1,1,1 };
			T weight_sum = 0, value_sum = 0;
			if constexpr (d == 2) {
				real w[2][2] = { {1.0 - frac[0],frac[0]},{1.0 - frac[1],frac[1]} };
				for (int s = 0; s < 4; s++) {
					int d0 = dx[s], d1 = dy[s];
					VectorDi refr_coord(coord[0] + d0, coord[1] + d1);
					if (grid.Valid(refr_coord) && valid_mask(refr_coord)) {
						T weight = w[0][d0] * w[1][d1];
						weight_sum += weight;
						value_sum += weight * data[grid.Index(coord[0] + d0, coord[1] + d1)];
					}
				}
			}
			else if constexpr (d == 3) {
				real w[3][2] = { {1.0 - frac[0],frac[0]},{1.0 - frac[1],frac[1]} ,{1.0 - frac[2],frac[2]} };
				for (int s = 0; s < 8; s++) {
					int d0 = dx[s], d1 = dy[s], d2 = dz[s];
					VectorDi refr_coord(coord[0] + d0, coord[1] + d1, coord[2] + d2);
					if (grid.Valid(refr_coord) && valid_mask(refr_coord)) {
						T weight = w[0][d0] * w[1][d1] * w[2][d2];
						weight_sum += weight;
						value_sum += weight * data[grid.Index(coord[0] + d0, coord[1] + d1, coord[2] + d2)];
					}
				}
			}
			else Assert("Interpolation:Linear_Intp error: dimension must be 2 or 3");
			return weight_sum == 0 ? 0 : value_sum / weight_sum;
		}
	};
}