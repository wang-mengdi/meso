#pragma once
#include "Common.h"
#include "AuxFunc.h"

namespace Meso {
	template<class T, int d, DataHolder side>
	T HausdorffDistance(const Array<Vector<T,d>, side>& a,const Array<Vector<T,d>, side>& b) {
		Typedef_VectorD(d);
		real sup_inf_a=0;
		//Can be parallelized
		for (int i = 0; i < b.size(); i++) {
			Array<VectorD,side> temp_a(a);
			VectorD b_i = b[i];
			VectorD n_b_i = - b_i;
			//Must write it this way because the add scalar function does not deal with the host array of vector well somehow
			if constexpr (side == HOST) {
				ArrayFunc::Add_Scalar(temp_a, -b_i);
			}
			else {
				ArrayFunc::Add_Scalar(temp_a, n_b_i);
			}
			real inf = ArrayFunc::Smallest_Norm<Vector<T,d>, side>(temp_a);
			if (inf > sup_inf_a) { sup_inf_a = inf; }
		}

		real sup_inf_b = 0;
		for (int i = 0; i < a.size(); i++) {
			Array<VectorD,side> temp_b(b);
			VectorD a_i = a[i];
			VectorD n_a_i = -a_i;
			if constexpr (side == HOST) {
				ArrayFunc::Add_Scalar(temp_b, -a_i);
			}
			else {
				ArrayFunc::Add_Scalar(temp_b, n_a_i);
			}
			real inf = ArrayFunc::Smallest_Norm<VectorD, side>(temp_b);
			if (inf > sup_inf_b) { sup_inf_b = inf; }
		}

		return std::max(sup_inf_a, sup_inf_b);
	}
}