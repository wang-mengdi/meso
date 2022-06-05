#pragma once
#include "SimplicialPrimitives.h"
#include "common.h"

using namespace Meso;
namespace ThinShellAuxFunc {
	inline int Element_Edges(const Vector2i& v, ArrayF<Vector<int, 1>, 2>& edges)
	{
		edges[0][0] = v[0]; edges[1][0] = v[1]; return 2;
	}

	inline int Element_Edges(const Vector3i& v, ArrayF<Vector<int, 2>, 3>& edges)
	{
		edges[0] = Vector2i(v[0], v[1]); edges[1] = Vector2i(v[1], v[2]); edges[2] = Vector2i(v[2], v[0]); return 3;
	}

	inline void Grad_Q(const ArrayF<Vector3, 3>& vtx, const int i, const int j, const Vector3& ps, const real& qs_i, const ArrayF<Vector3, 3>& ls, const real& a, Vector3& grad_q) {
		int jr = (j + 1) % 3;
		int jl = (j + 2) % 3;

		real tmp1 = ps[jr] * qs_i / (real)8 / a / a;
		if (i == jr) { tmp1 -= (real)0.5 / (a * ls[i].norm()); }

		real tmp2 = ps[jl] * qs_i / (real)8 / a / a;
		if (i == jl) { tmp2 -= (real)0.5 / (a * ls[i].norm()); }

		grad_q = tmp1 * ls[jr] - tmp2 * ls[jl];
	}

	inline void Grad_R(const ArrayF<Vector3, 3>& vtx, const int i, const int j, const Vector3& ps, const ArrayF<Vector3, 3>& ls, const real& a, ArrayF<Vector3, 2>& grad_r) {
		int jr = (j + 1) % 3;
		int jl = (j + 2) % 3;

		Vector3 grad_a = (real)1 / (real)8 / a * (ps[jl] * ls[jl] - ps[jr] * ls[jr]);
		Vector3 grad_l = ((int)(i == jl)) / ls[jl].norm() * ls[jl] - ((int)(i == jr)) / ls[jr].norm() * ls[jr];

		Vector3 coef2 = (grad_a * ls[i].norm() + grad_l * a) / ((real)4 * a * a * ls[i].dot(ls[i]));
		real coef1 = a * ls[i].norm() * (real)2;

		int i_p = (i + 2) % 3;
		Vector3 tmp1 = ((((int)(i_p == jr)) + ((int)(i_p == j)) - ((int)(i_p == jl))) * ls[jl] + (((int)(i_p == jr)) - ((int)(i_p == j)) - ((int)(i_p == jl))) * ls[jr]);
		grad_r[0] = tmp1 / coef1 - ps[i_p] * coef2;

		i_p = (i + 1) % 3;
		tmp1 = ((((int)(i_p == jr)) + ((int)(i_p == j)) - ((int)(i_p == jl))) * ls[jl] + (((int)(i_p == jr)) - ((int)(i_p == j)) - ((int)(i_p == jl))) * ls[jr]);
		grad_r[1] = tmp1 / coef1 - ps[i_p] * coef2;
	}

	inline void Grad_N(const ArrayF<Vector3, 3>& vtx, const ArrayF<Vector3, 3>& ls, const real& a, ArrayF<Matrix3, 3>& grad_n) {
		Vector3 n = Triangle<3>::Normal(vtx[0], vtx[1], vtx[2]);
		for (int i = 0; i < 3; i++) {
			grad_n[i] = -n * ((real)0.5 * ls[i].cross(n) / a).transpose();
		}
	}
};
