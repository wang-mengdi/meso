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
			if (iter_num == 0) {
				T* x_ptr = thrust::raw_pointer_cast(x.data());
				cudaMemset(x_ptr, 0, sizeof(T) * dof);
				return;
			}
			{
				T* x_temp_ptr = thrust::raw_pointer_cast(x_temp.data());
				cudaMemset(x_temp_ptr, 0, sizeof(T) * dof);
			}
			for (int i = 0; i < iter_num; i++) {
				//b-Ax
				mapping->Residual(x, x_temp, b);
				//(b-Ax)/.diag
				T* x_ptr = thrust::raw_pointer_cast(x.data());
				T* one_over_diag_ptr = thrust::raw_pointer_cast(one_over_diag.data());
				auto mul = [=]__device__(T & a, T & b) { a *= b; };
				GPUFunc::Cwise_Mapping_Wrapper(x_ptr, one_over_diag_ptr, mul, dof);
				//x+=(b-Ax)/.diag*.omega
				real _omega = omega;
				T* x_temp_ptr = thrust::raw_pointer_cast(x_temp.data());
				auto weighted_add = [=]__device__(T & a, T & b) { a = b + a * _omega; };
				GPUFunc::Cwise_Mapping_Wrapper(x_ptr, x_temp_ptr, weighted_add, dof);
				if (i + 1 < iter_num)cudaMemcpy(x_temp_ptr, x_ptr, sizeof(T) * dof, cudaMemcpyDeviceToDevice);
			}
			checkCudaErrors(cudaGetLastError());
		}
	};
}
