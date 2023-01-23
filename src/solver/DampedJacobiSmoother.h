//////////////////////////////////////////////////////////////////////////
// Damped Jacobian Smoother
// Copyright (c) (2022-), Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once

#include "PoissonMapping.h"

namespace Meso {
	template<class T>
	__global__ void Jacobi_Boundary_Smooth(const int* _boundary_tiles, const unsigned char* _cell_type, const T* _b, 
		const T* _Ax, const T* _one_over_diag, const T _omega, T* _x)
	{
		int cell_id = _boundary_tiles[blockIdx.x] * 64 + threadIdx.x;
		if (_cell_type[cell_id] != 3)
			return;
		_x[cell_id] += (_b[cell_id] - _Ax[cell_id]) * _one_over_diag[cell_id] * _omega;
	}

	//with zero initial guess, so it's a linear mapping of b in Ax=b
	template<class T, int d>
	class DampedJacobiSmoother : public LinearMapping<T> {
	public:
		using Base=LinearMapping<T>;
		MaskedPoissonMapping<T, d>* poisson;
		T omega;
		int dof;
		int iter_num;
		int boundary_iter_num;
		ArrayDv<T> one_over_diag;
		ArrayDv<T> x_temp;
		DampedJacobiSmoother() {}
		DampedJacobiSmoother(MaskedPoissonMapping<T, d>& _poisson, const ArrayDv<T>& _diag, const int _iter_num, 
			const int _boundary_iter_num, const real _omega = 2.0 / 3.0) { Init(_poisson, _diag, _iter_num, _boundary_iter_num, _omega); }
		void Init(MaskedPoissonMapping<T, d>& _poisson, const ArrayDv<T>& _diag, const int _iter_num, 
			const int _boundary_iter_num, const T _omega = 2.0 / 3.0) {
			poisson = &_poisson;
			iter_num = _iter_num;
			boundary_iter_num = _boundary_iter_num;
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
		void Boundary_Apply(ArrayDv<T>& x, const ArrayDv<T>& b)
		{
			for (int i = 0; i < boundary_iter_num; i++)
			{
				//Ax
				poisson->Boundary_Apply(x_temp, x);
				//x+=b-Ax/.diag*.omega
				const int* _boundary_tiles; const unsigned char* _cell_type; const T* _b;
				const T* _Ax; const T* _one_over_diag; const T _omega = 0; T* _x;
				Jacobi_Boundary_Smooth<T> << < poisson->boundary_tiles.size(), 64 >> > (ArrayFunc::Data(poisson->boundary_tiles), 
					poisson->cell_type.Data_Ptr(), ArrayFunc::Data(b),ArrayFunc::Data(x_temp), 
					ArrayFunc::Data(one_over_diag), omega, ArrayFunc::Data(x));
			}
		}
	};
}
