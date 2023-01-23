//////////////////////////////////////////////////////////////////////////
// Damped Jacobian Smoother
// Copyright (c) (2022-), Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once

#include "PoissonMapping.h"

namespace Meso {
	//with zero initial guess, so it's a linear mapping of b in Ax=b
	template<class T, int d>
	class DampedJacobiSmoother : public LinearMapping<T> {
	public:
		using Base=LinearMapping<T>;
		MaskedPoissonMapping<T, d>* poisson;
		T omega;
		int dof;
		int iter_num;
		ArrayDv<T> one_over_diag;
		ArrayDv<T> x_temp;
		DampedJacobiSmoother() {}
		DampedJacobiSmoother(MaskedPoissonMapping<T, d>& _poisson, const ArrayDv<T>& _diag, const int _iter_num, const real _omega = 2.0 / 3.0) { Init(_poisson, _diag, _iter_num, _omega); }
		void Init(MaskedPoissonMapping<T, d>& _poisson, const ArrayDv<T>& _diag, const int _iter_num, const T _omega = 2.0 / 3.0) {
			poisson = &_poisson;
			iter_num = _iter_num;
			omega = _omega;
			dof = poisson->XDoF();
			one_over_diag.resize(dof);
			const T* diag_ptr = thrust::raw_pointer_cast(_diag.data());
			T* one_over_diag_ptr = thrust::raw_pointer_cast(one_over_diag.data());
			auto divide = [=]__device__(T & a, const T & b) { a = 1.0 / b; };
			GPUFunc::Cwise_Mapping_Wrapper(one_over_diag_ptr, diag_ptr, divide, dof);
			x_temp.resize(dof);
		}
		virtual int XDoF()const { return dof; }
		virtual int YDoF()const { return dof; }
		virtual void Apply(ArrayDv<T>& x, const ArrayDv<T>& b) {
			Base::Memory_Check(x, b, "DampedJacobiSmoother::Apply error: not enough memory space");
			for (int i = 0; i < iter_num; i++) {
				//b-Ax
				poisson->Residual(x_temp, x, b);
				//x+=(b-Ax)/.diag*.omega
				T _omega = omega;
				auto func = [_omega]__device__(T & a, const T & b, const T & c) { a += b * c * _omega; };
				GPUFunc::Cwise_Mapping_Wrapper(ArrayFunc::Data(x), ArrayFunc::Data(x_temp), 
					ArrayFunc::Data(one_over_diag), func, dof);
			}
			checkCudaErrors(cudaGetLastError());
		}
	};
}
