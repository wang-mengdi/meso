#pragma once
#include "Common.h"
#include "AuxFunc.h"

namespace Meso {
	template<class T, int d>
	T HausdorffDistance(const Array<Vector<T,d>>& a,const Array<Vector<T,d>>& b) {
		real sup_inf_a=0;
		for (int i = 0; i < b.size(); i++) {
			Array<Vector<T, d>> temp_a = a;
			ArrayFunc::Add_Scalar(temp_a, -b[i]);			
			real inf = ArrayFunc::Smallest_Norm<Vector<T,d>>(temp_a);
			if (inf > sup_inf_a) { sup_inf_a = inf; }
		}

		real sup_inf_b = 0;
		for (int i = 0; i < a.size(); i++) {
			Array<Vector<T, d>> temp_b = b;
			ArrayFunc::Add_Scalar(temp_b, -a[i]);
			real inf = ArrayFunc::Smallest_Norm<Vector<T, d>>(temp_b);
			if (inf > sup_inf_b) { sup_inf_b = inf; }
		}

		return std::max(sup_inf_a, sup_inf_b);
	}
}