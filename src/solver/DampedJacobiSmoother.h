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
				//(b-Ax)/.diag
				T* x_temp_ptr = ArrayFunc::Data(x_temp);
				T* one_over_diag_ptr = ArrayFunc::Data(one_over_diag);
				auto mul = [=]__device__(T & a, T & b) { a *= b; };
				GPUFunc::Cwise_Mapping_Wrapper(x_temp_ptr, one_over_diag_ptr, mul, dof);
				//x+=(b-Ax)/.diag*.omega
				real _omega = omega;
				T* x_ptr = ArrayFunc::Data(x);
				auto weighted_add = [=]__device__(T & a, T & b) { a += b * _omega; };
				GPUFunc::Cwise_Mapping_Wrapper(x_ptr, x_temp_ptr, weighted_add, dof);
			}
			checkCudaErrors(cudaGetLastError());
		}
	};
}
