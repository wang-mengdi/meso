//////////////////////////////////////////////////////////////////////////
// Damped Jacobian Smoother
// Copyright (c) (2022-), Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once

#include "LinearMapping.h"
#include "PoissonFunc.h"

namespace Meso {
	//with zero initial guess, so it's a linear mapping of b in Ax=b
	template<class T>
	class DampedJacobiSmoother : public LinearMapping<T> {
	public:
		LinearMapping<T>* mapping;
		T omega;
		int dof;
		int iter_num;
		ArrayDv<T> diag;
		ArrayDv<T> x_temp;
		DampedJacobiSmoother() {}
		template<int d> DampedJacobiSmoother(MaskedPoissonMapping<T, d>& _mapping, const int _iter_num, const real _omega = 2.0 / 3.0) { Init(_mapping, _iter_num, _omega); }
		template<int d>
		void Init(MaskedPoissonMapping<T, d>& _mapping, const int _iter_num, const T _omega = 2.0 / 3.0) {
			mapping = &_mapping;
			iter_num = _iter_num;
			omega = _omega;
			dof = mapping->XDof();
			Poisson_Diagonal(diag, _mapping);
			x_temp.resize(dof);
		}
		virtual int XDof()const { return dof; }
		virtual int YDof()const { return dof; }
		virtual void Apply(ArrayDv<T>& x, const ArrayDv<T>& b) {
			Memory_Check(x, b, "DampedJacobiSmoother::Apply error: not enough memory space");
			if (iter_num == 0) {
				ArrayFunc::Fill(x, (T)0);
				return;
			}
			ArrayFunc::Fill(x_temp, (T)0);
			for (int i = 0; i < iter_num; i++) {
				//b-Ax
				mapping->Residual(x, x_temp, b);
				//(b-Ax)/.diag
				ArrayFunc::Binary_Transform(x, diag, [=]__device__(T a, T b) { return a / b; }, x);
				//x+=(b-Ax)/.diag*.omega
				real _omega = omega;
				ArrayFunc::Binary_Transform(x, x_temp, [=]__device__(T a, T b) { return b + a * _omega; }, x);

				if (i + 1 < iter_num)ArrayFunc::Copy(x_temp, x);
			}
			checkCudaErrors(cudaGetLastError());
		}
	};
}