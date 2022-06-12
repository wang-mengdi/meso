//////////////////////////////////////////////////////////////////////////
// Damped Jacobian Smoother
// Copyright (c) (2022-), Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once

#include "LinearMapping.h"

namespace Meso {
	//with zero initial guess, so it's a linear mapping of b in Ax=b
	template<class T>
	class DampedJacobiSmoother : public LinearMapping<T> {
	public:
		using Base=LinearMapping<T>;
		LinearMapping<T>* mapping;
		T omega;
		int dof;
		int iter_num;
		ArrayDv<T> one_over_diag;
		ArrayDv<T> x_temp;
		DampedJacobiSmoother() {}
		DampedJacobiSmoother(LinearMapping<T>& _mapping, const ArrayDv<T>& _diag, const int _iter_num, const real _omega = 2.0 / 3.0) { Init(_mapping, _diag, _iter_num, _omega); }
		void Init(LinearMapping<T>& _mapping, const ArrayDv<T>& _diag, const int _iter_num, const T _omega = 2.0 / 3.0) {
			mapping = &_mapping;
			iter_num = _iter_num;
			omega = _omega;
			dof = mapping->XDoF();
			one_over_diag.resize(dof);
			ArrayFunc::Unary_Transform(_diag, 1.0 / thrust::placeholders::_1, one_over_diag);
			x_temp.resize(dof);
		}
		virtual int XDoF()const { return dof; }
		virtual int YDoF()const { return dof; }
		virtual void Apply(ArrayDv<T>& x, const ArrayDv<T>& b) {
			Base::Memory_Check(x, b, "DampedJacobiSmoother::Apply error: not enough memory space");
			if (iter_num == 0) {
				ArrayFunc::Fill(x, (T)0);
				return;
			}
			ArrayFunc::Fill(x_temp, (T)0);
			for (int i = 0; i < iter_num; i++) {
				//b-Ax
				mapping->Residual(x, x_temp, b);
				//(b-Ax)/.diag
				ArrayFunc::Multiply(x, one_over_diag);
				//x+=(b-Ax)/.diag*.omega
				real _omega = omega;
				ArrayFunc::Binary_Transform(x, x_temp, [=]__device__(T a, T b) { return b + a * _omega; }, x);

				if (i + 1 < iter_num)ArrayFunc::Copy(x_temp, x);
			}
			checkCudaErrors(cudaGetLastError());
		}
	};
}
