//////////////////////////////////////////////////////////////////////////
// Rendering functions
// Copyright (c) (2022-), Yitong Deng
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "Kernels.h"
#include "ParticleFrames.h"

namespace Meso {

	template<int d, class T>
	T Surface_Sum_SPH(const Vector<real,d>& my_pos, const Matrix<real, d>& frame, 
		const Array<T>& f, const Array<Vector<real,d>>& pos, const Array<int>& nbs,
		const KernelSPH& kernel, const real radius, KernelType kernel_type = KernelType::QUINTIC) {
		Typedef_VectorD(d);
		T sum = VectorFunc::Zero<T>();
		for (int k = 0; k < nbs.size(); k++) {
			int j = nbs[k];
			VectorD wr = my_pos - pos[j];
			VectorT sr = Rotate_To_TPlane<d>(wr, frame);
			real w = kernel.Weight<d - 1>(sr, radius, kernel_type);
			sum += f[j] * w;
		}
		return sum;
	}

	template<int d, class T>
	T Surface_Laplacian_SPH(const Vector<real, d>& my_pos, const Matrix<real, d>& frame,
		const std::function<T(const int)>& f, const T& my_f, const Array<real>& a,
		const Array<Vector<real, d>>& pos, const Array<int>& nbs,
		const KernelSPH& kernel, const real radius, KernelType kernel_type = KernelType::QUINTIC) {
		Typedef_VectorD(d);
		T lap = VectorFunc::Zero<T>();
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

	template<int d>
	Vector<real,d-1> Surface_Gradient_Diff_SPH(const Vector<real, d>& my_pos, const Matrix<real, d>& frame,
		const Array<real>& f, const real& my_f, const Array<real>& a,
		const Array<Vector<real, d>>& pos, const Array<int>& nbs,
		const KernelSPH& kernel, const real radius, KernelType kernel_type = KernelType::QUINTIC) {
		Typedef_VectorD(d);
		VectorT grad_f = VectorT::Zero();
		for (int k = 0; k < nbs.size(); k++) {
			int j = nbs[k];
			VectorD wr = my_pos - pos[j];
			VectorT sr = Rotate_To_TPlane<d>(wr, frame);
			VectorT grad_w = kernel.Grad<d - 1>(sr, radius, kernel_type);
			grad_f += a[j] * (f[j] - my_f) * grad_w;
		}
		return grad_f;
	}
}