//////////////////////////////////////////////////////////////////////////
// Particle Frames
// Copyright (c) (2018-), Bo Zhu, Hui Wang, Mengdi Wang, Yitong Deng
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once

#include "Common.h"
#include "AuxFunc.h"
#include "IOFunc.h"

namespace Meso {

	////2D->3D projection, back to world space
	inline void Unproject_To_World(const Vector2& t, const Matrix3& e, Vector3& u) { u = e.col(0) * t[0] + e.col(1) * t[1]; }
	////1D->2D projection, back to world space
	inline void Unproject_To_World(const Vector1& t, const Matrix2& e, Vector2& u) { u = e.col(0) * t[0]; }

	template<int d>
	static Vector<real, d> Project_To_Norm(const Vector<real, d>& u, const Matrix<real, d>& E) {
		return u.dot(E.col(d - 1)) * E.col(d - 1);
	}

	template<int d>
	static Vector<real, d - 1> Project_To_TPlane(const Vector<real, d>& u, const Matrix<real, d>& E) {
		Vector<real, d - 1> t_coords;
		for (int i = 0; i < d - 1; i++) {
			t_coords[i] = u.dot(E.col(i));
		}
		return t_coords;
	}
	template<int d>
	static Vector<real, d - 1> Rotate_To_TPlane(const Vector<real, d>& u, const Matrix<real, d>& E) {//same as Project_To_TPlane, but preserves length
		Vector<real, d - 1> t_coords;
		for (int i = 0; i < d - 1; i++) {
			t_coords[i] = u.dot(E.col(i));
		}
		if (t_coords.norm() == 0.) return Vector<real, d - 1>::Zero();
		else return u.norm() * t_coords.normalized();
	}

	////Kernel for PCA normal calculation
	inline real W_PCA(const real d, const real r)
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
		VectorD normal = MathFunc::Min_Eigenvector(C);

		// Always align with previous normal
		if (normal.dot(my_frame.col(d-1)) < (real) 0) normal *= (real)-1;

		////update local frame according to the PCA normal
		if constexpr (d == 2) {
			VectorD tang = -MathFunc::Orthogonal_Vector(normal);
			my_frame.col(0) = tang.normalized();
			my_frame.col(1) = normal.normalized();
		}
		else if constexpr (d == 3) {
			VectorD t1 = -MathFunc::Orthogonal_Vector(normal);
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