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
	class SPH {
		Typedef_VectorD(d); Typedef_MatrixD(d);
	private:
		KernelSPH kernel;

	public:
		SPH() {}	

		template<class T>
		T Sum(const VectorD& my_pos, const MatrixD& frame,
			const std::function<T(const int)>& f, const Array<VectorD>& pos, const Array<int>& nbs,
			const real radius, KernelType kernel_type = KernelType::QUINTIC) {
			Typedef_VectorD(d);
			T sum = MathFunc::Zero<T>();
			for (int k = 0; k < nbs.size(); k++) {
				int j = nbs[k];
				VectorD wr = my_pos - pos[j];
				VectorT sr = Rotate_To_TPlane<d>(wr, frame);
				real w = kernel.Weight<d - 1>(sr, radius, kernel_type);
				sum += f(j) * w;
			}
			return sum;
		}

		template<class T>
		T Sum(const VectorD& my_pos, const MatrixD& frame,
			const Array<T>& f, const Array<VectorD>& pos, const Array<int>& nbs,
			const real radius, KernelType kernel_type = KernelType::QUINTIC) {
			std::function<T(const int)>& f_func = [&](const int idx)->T {return f[idx];};
			return Sum(my_pos, frame, f_func, pos, nbs, radius, kernel_type);
		}

		template<class T>
		T Laplacian(const VectorD& my_pos, const MatrixD& frame,
			const std::function<T(const int)>& f, const T& my_f, const Array<real>& a,
			const Array<VectorD>& pos, const Array<int>& nbs,
			const real radius, KernelType kernel_type = KernelType::QUINTIC) {
			T lap = MathFunc::Zero<T>();
			for (int k = 0; k < nbs.size(); k++) {
				int j = nbs[k];
				VectorD wr = my_pos - pos[j];
				VectorT sr = Rotate_To_TPlane<d>(wr, frame);
				real norm = std::max(sr.norm(), radius * 0.001);
				VectorT grad_w = kernel.Grad<d - 1>(sr, radius, kernel_type);
				lap += a[j] * (f(j) - my_f) * 2 * grad_w.norm() / norm;
			}
			return lap;
		}

		template<class T>
		T Laplacian(const VectorD& my_pos, const MatrixD& frame,
			const Array<T>& f, const T& my_f, const Array<real>& a,
			const Array<VectorD>& pos, const Array<int>& nbs,
			const real radius, KernelType kernel_type = KernelType::QUINTIC) {
			std::function<T(const int)>& f_func = [&](const int idx)->T {return f[idx]; };
			return Laplacian(my_pos, frame, f_func, my_f, a, pos, nbs, radius, kernel_type);
		}

		real Divergence(const VectorD& my_pos, const MatrixD& frame,
			const std::function<VectorT(const int)>& f, const VectorT& my_f, const Array<real>& a,
			const Array<VectorD>& pos, const Array<int>& nbs,
			const real radius, KernelType kernel_type = KernelType::QUINTIC) {
			real div = 0.;
			for (int k = 0; k < nbs.size(); k++) {
				int j = nbs[k];
				VectorD wr = my_pos - pos[j];
				VectorT sr = Rotate_To_TPlane<d>(wr, frame);
				VectorT grad_w = kernel.Grad<d - 1>(sr, radius, kernel_type);
				VectorT fr = f(j) - my_f;
				div += a[j] * fr.dot(grad_w);
			}
			return div;
		}

		real Divergence(const VectorD& my_pos, const MatrixD& frame,
			const Array<VectorT>& f, const VectorT& my_f, const Array<real>& a,
			const Array<VectorD>& pos, const Array<int>& nbs,
			const real radius, KernelType kernel_type = KernelType::QUINTIC) {
			std::function<VectorT(const int)>& f_func = [&](const int idx)->VectorT {return f[idx]; };
			return Divergence(my_pos, frame, f_func, my_f, a, pos, nbs, radius, kernel_type);
		}

		VectorT Gradient_Diff(const VectorD& my_pos, const MatrixD& frame,
			const std::function<real(const int)>& f, const real& my_f, const Array<real>& a,
			const Array<VectorD>& pos, const Array<int>& nbs,
			const real radius, KernelType kernel_type = KernelType::QUINTIC) {
			VectorT grad_f = VectorT::Zero();
			for (int k = 0; k < nbs.size(); k++) {
				int j = nbs[k];
				VectorD wr = my_pos - pos[j];
				VectorT sr = Rotate_To_TPlane<d>(wr, frame);
				VectorT grad_w = kernel.Grad<d - 1>(sr, radius, kernel_type);
				grad_f += a[j] * (f(j) - my_f) * grad_w;
			}
			return grad_f;
		}

		VectorT Gradient_Diff(const VectorD& my_pos, const MatrixD& frame,
			const Array<real>& f, const real& my_f, const Array<real>& a,
			const Array<VectorD>& pos, const Array<int>& nbs,
			const real radius, KernelType kernel_type = KernelType::QUINTIC) {
			std::function<real(const int)>& f_func = [&](const int idx)->real {return f[idx]; };
			return Gradient_Diff(my_pos, frame, f_func, my_f, a, pos, nbs, radius, kernel_type);
		}
	};

}