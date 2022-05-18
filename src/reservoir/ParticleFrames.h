//////////////////////////////////////////////////////////////////////////
// Particle Frames
// Copyright (c) (2018-), Bo Zhu, Hui Wang, Mengdi Wang, Yitong Deng
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////

#include "Common.h"
#include "AuxFunc.h"
#include "IOFunc.h"

namespace Meso {

	////Kernel for PCA normal calculation
	real W_PCA(const real d, const real r)
	{
		if (d < r)return (real)1 - pow(d / r, 3);
		else return (real)0;
	}

	template<int d>
	void Set_Local_Frame_PCA(const real r, const Vector<real, d>& my_pos, const Array<Vector<real, d>>& nbs_pos, Matrix<real, d>& my_frame) {
		Typedef_VectorD(d); Typedef_MatrixD(d);
		size_t nbs_num = nbs_pos.size();
		VectorD xp = VectorD::Zero();
		real w = (real)0;
		for (size_t q = 0; q < nbs_num; q++) {
			////here we use volumetric distance instead of tangential distance
			real dis = (nbs_pos[q] - my_pos).norm();
			real w0 = W_PCA(dis, r);
			xp += w0 * nbs_pos[q];
			w += w0;
		}
		if (w != (real)0)xp /= w;
		MatrixD C = MatrixD::Zero();
		real wc = (real)0;
		for (size_t q = 0; q < nbs_num; q++) {
			const VectorD& xq = nbs_pos[q];
			real dis = (xq - xp).norm();
			real w0 = W_PCA(dis, r);
			C += w0 * (xq - xp) * (xq - xp).transpose();
			wc += w0;
		}
		if (wc != (real)0)C /= wc;
		VectorD normal = VectorFunc::Min_Eigenvector(C);

		// Always align with previous normal
		if (normal.dot(my_frame.col(d-1)) < (real) 0) normal *= (real)-1;

		////update local frame according to the PCA normal
		if constexpr (d == 2) {
			VectorD tang = -VectorFunc::Orthogonal_Vector(normal);
			my_frame.col(0) = tang.normalized();
			my_frame.col(1) = normal.normalized();
		}
		else if constexpr (d == 3) {
			VectorD t1 = -VectorFunc::Orthogonal_Vector(normal);
			VectorD t2 = t1.cross(normal);
			my_frame.col(0) = t1.normalized();
			my_frame.col(1) = t2.normalized();
			my_frame.col(2) = normal.normalized();
		}
	}

	// Orient normal directions
	// Pointing Outward from COM, can be used on simple shapes like box, or sphere.
	template<int d>
	void Orient_Normals_COM(const Array<Vector<real, d>>& poss, Array<Matrix<real, d>>& frames) {
		Typedef_VectorD(d); Typedef_MatrixD(d);
		Assert(poss.size() == frames.size(), "[Partialce Frames] positions and frames have different sizes {}, vs {}", poss.size(), frames.size());
		VectorD COM = ArrayFunc::Mean<VectorD>(poss);
#pragma omp parallel for
		for (int i = 0; i < frames.size(); i++) {
			VectorD outward = poss[i] - COM;
			MatrixD& my_frame = frames[i];
			if (outward.dot(my_frame.col(d - 1)) < (real)0) {
				my_frame.col(d - 1) *= (real)-1;
				if (d > 1) my_frame.col(0) *= (real)-1; //don't forget to flip another axis (not considering d > 3)
			}
		}
	}
}