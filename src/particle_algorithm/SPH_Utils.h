//////////////////////////////////////////////////////////////////////////
// SPH
// Copyright (c) (2022-), Yitong Deng
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "Kernels.h"
#include "ParticleFrames.h"

namespace Meso {

	template<int d>
	class SPH_Utils {
		Typedef_VectorD(d); Typedef_MatrixD(d);
	public:
		KernelSPH kernel;

	public:
		SPH_Utils() {}	

		template<class IFuncT, class IFuncV>
		decltype(auto) Sum(IFuncT val_func, IFuncV diff_func, const Array<int>& nbs,
			const real radius, KernelType kernel_type = KernelType::QUINTIC) {
			using T = decltype(val_func(0));
			using V = decltype(diff_func(0));
			T sum = MathFunc::Zero<T>();
			for (int k = 0; k < nbs.size(); k++) {
				int j = nbs[k];
				V wr = diff_func(j);
				real w = kernel.Weight<d>(wr, radius, kernel_type);
				sum += val_func(j) * w;
			}
			return sum;
		}

		template<class IFuncT, class IFuncV>
		decltype(auto) Sum_Grad(IFuncT val_func, IFuncV diff_func, const Array<int>& nbs,
			const real radius, KernelType kernel_type = KernelType::QUINTIC) {
			using T = decltype(val_func(0));
			using V = decltype(diff_func(0));
			VectorD sum_grad = VectorD::Zero();
			for (int k = 0; k < nbs.size(); k++) {
				int j = nbs[k];
				V wr = diff_func(j);
				VectorD grad_w = kernel.Grad<d>(wr, radius, kernel_type);
				sum_grad += val_func(j) * grad_w;
			}
			return sum_grad;
		}

		template<class IFuncT, class IFuncV, class IFuncA>
		decltype(auto) Grad(IFuncT val_func, IFuncV diff_func, IFuncA vol_func, 
			const Array<int>& nbs, const real radius, KernelType kernel_type = KernelType::QUINTIC) {
			using T = decltype(val_func(0));
			using V = decltype(diff_func(0));
			V grad_f = MathFunc::Zero<V>();
			for (int k = 0; k < nbs.size(); k++) {
				int j = nbs[k];
				V wr = diff_func(j);
				VectorD grad_w = kernel.Grad<d>(wr, radius, kernel_type);
				grad_f += vol_func(j) * val_func(j) * grad_w;
			}
			return grad_f;
		}

		template<class IFuncT, class IFuncV, class IFuncA>
		decltype(auto) Lap_Brookshaw (IFuncT val_func, IFuncV diff_func, IFuncA vol_func,
			const Array<int>& nbs, const real radius, KernelType kernel_type = KernelType::QUINTIC) {
			using T = decltype(val_func(0));
			using V = decltype(diff_func(0));
			T lap = MathFunc::Zero<T>();
			for (int k = 0; k < nbs.size(); k++) {
				int j = nbs[k];
				V wr = diff_func(j);
				real norm = std::max(wr.norm(), radius * 0.001);
				VectorD grad_w = kernel.Grad<d>(wr, radius, kernel_type);
				lap += vol_func(j) * val_func(j) * 2 * grad_w.norm() / norm;
			}
			return lap;
		}

		template<class IFuncV, class IFuncA>
		decltype(auto) Lap_Diagonal_Diff(const real vol_i, IFuncV diff_func, IFuncA vol_func,
			const Array<int>& nbs, const real radius, KernelType kernel_type = KernelType::QUINTIC) {
			using V = decltype(diff_func(0));
			using T = decltype(vol_func(0));
			VectorD temp = VectorD::Zero();
			for (int k = 0; k < nbs.size(); k++) {
				int j = nbs[k];
				V wr = diff_func(j);
				VectorD grad_w = kernel.Grad<d>(wr, radius, kernel_type);
				temp += vol_func(j) * grad_w;
			}
			real a_ii = 0.;
			for (int k = 0; k < nbs.size(); k++) {
				int j = nbs[k];
				V wr = diff_func(j);
				VectorD grad_w = kernel.Grad<d>(wr, radius, kernel_type);
				a_ii += vol_func(j) * (-vol_i * grad_w + 0. * temp).dot(grad_w);
			}
			return a_ii;
		}

		template<class IFuncT, class IFuncV, class IFuncA>
		decltype(auto) Div(IFuncT val_func, IFuncV diff_func, IFuncA vol_func,
			const Array<int>& nbs, const real radius, KernelType kernel_type = KernelType::QUINTIC) {
			using V = decltype(diff_func(0));
			using T = V; // in div computation, the summand must be in the same dimension as the position difference
			real div = 0.;
			for (int k = 0; k < nbs.size(); k++) {
				int j = nbs[k];
				V wr = diff_func(j);
				VectorD grad_w = kernel.Grad<d>(wr, radius, kernel_type);
				div += vol_func(j) * val_func(j).dot(grad_w);
			}
			return div;
		}
	};

}