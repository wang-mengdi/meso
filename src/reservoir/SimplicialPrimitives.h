#pragma once
#include "Common.h"
#include "AuxFunc.h"

namespace Meso{
	//////////////////////////////////////////////////////////////////////////
	////Segment
	//////////////////////////////////////////////////////////////////////////

	template<int d> class Segment
	{
		Typedef_VectorD(d);
	public:
		VectorD p0;
		VectorD p1;

		Segment() :p0(), p1() {}
		Segment(const VectorD& _p0, const VectorD& _p1) :p0(_p0), p1(_p1) {}
	};


	//////////////////////////////////////////////////////////////////////////
	////Triangle
	//////////////////////////////////////////////////////////////////////////

	template<int d> class Triangle{
		Typedef_VectorD(d);
	public:
		VectorD p0, p1, p2;
		Triangle(const VectorD& _p0, const VectorD& _p1, const VectorD& _p2) :p0(_p0), p1(_p1), p2(_p2) {}

		virtual bool Inside(const VectorD& pos) const { 
			if constexpr (d == 2) { return Inside(p0, p1, p2, pos); }
			else if constexpr (d == 3) { return false; }
		}

		virtual VectorD Normal(const VectorD& pos) const { 
			if constexpr (d == 3) { return Normal(p0, p1, p2); }
			else { Error("Dimension not supported"); }
		}

		static bool Inside(const VectorD& p0, const VectorD& p1, const VectorD& p2, const VectorD& pos)
		{
			if constexpr (d == 2) {
				VectorD v0 = p1 - p0; VectorD q0 = pos - p0; real c0 = MathFunc::Cross(v0, q0)[0];
				VectorD v1 = p2 - p1; VectorD q1 = pos - p1; real c1 = MathFunc::Cross(v1, q1)[0];
				VectorD v2 = p0 - p2; VectorD q2 = pos - p2; real c2 = MathFunc::Cross(v2, q2)[0];
				return (c0 > (real)0 && c1 > (real)0 && c2 > (real)0) || (c0 < (real)0 && c1 < (real)0 && c2 < (real)0);
			}
			else {
				Error("Dimension not supported");
			}
		}
	
		static int Third_Vertex(int p0, int p1, const Vector3i& vtxs) {
			for (int i = 0; i < 3; i++) {
				if (p0 != vtxs[i] && p1 != vtxs[i]) { return vtxs[i]; }
			}
			Error("Invalid input of existing indices");
			return -1; //shouldn't reach here
		}

		static VectorD Normal(const VectorD& p0, const VectorD& p1, const VectorD& p2)
		{
			if constexpr(d == 3){
				VectorD e0 = p1 - p0; VectorD e1 = p2 - p0; return e0.cross(e1).normalized();
			}
			else {
				Error("Dimension not supported");
			}
		}

		static real Area(const VectorD& p0, const VectorD& p1, const VectorD& p2)
		{
			if constexpr (d == 3) {
				Vector3 e0 = p1 - p0; Vector3 e1 = p2 - p0;
				return (real)0.5 *e0.cross(e1).norm();
			}
			else {
				Error("Dimension not supported");
			}
		}
	};
}