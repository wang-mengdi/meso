#pragma once
#include "AuxFunc.h"

namespace Meso {
	namespace WaveFunc {

		C Vector2C_Dot(const Vector2C& a, const Vector2C& b);

		Vector4 Vector2C_To_Vector4(const Vector2C& v);

		Vector2C Vector4_To_Vector2C(const Vector4& v);

		real Psi_To_Vel(Vector4 p0, Vector4 p1, const real h_bar, const real dx);

		template<int d> 
		Vector2C Vel_To_Psi_C(const Vector<real, d>& vel, const Vector<real, d>& pos) {
			Vector2C psi; psi[0] = C(1., 0.); psi[1] = C(.1, 0.);
			real norm = sqrt(thrust::norm(psi[0]) + thrust::norm(psi[1]));
			psi[0] /= norm; psi[1] /= norm;
			real phase = vel.dot(pos);
			for (int i = 0; i < 2; i++) { psi[i] *= exp(i * phase); }
			return psi;
		}
	}

}