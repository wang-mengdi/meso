#include "WaveFunc.h"

namespace Meso {
	namespace WaveFunc {
		C Vector2C_Dot(const Vector2C& a, const Vector2C& b)
		{
			return a[0] * b[0] + a[1] * b[1];
		}

		Vector4 Vector2C_To_Vector4(const Vector2C& v)
		{
			return Vector4(v[0].real(), v[0].imag(), v[1].real(), v[1].imag());
		}

		Vector2C Vector4_To_Vector2C(const Vector4& v)
		{
			Vector2C c(C(v[0], v[1]), C(v[2], v[3])); return c;
		}
		
		real Psi_To_Vel(Vector4 p0, Vector4 p1, const real h_bar, const real dx)
		{
			p0[1] *= -1;
			p0[3] *= -1;

			Vector2C psi_0 = Vector4_To_Vector2C(p0);
			Vector2C psi_1 = Vector4_To_Vector2C(p1);

			C q = psi_0.dot(psi_1);
			real psi_vel = thrust::arg(q) * h_bar / dx;
			return psi_vel;
		}
	}
}