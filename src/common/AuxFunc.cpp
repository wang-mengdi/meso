#include "AuxFunc.h"

namespace Meso {

	namespace VectorFunc {
		template<> Vector1 V<1>(const real x, const real y, const real z) { return Vector1(x); }
		template<> Vector2 V<2>(const real x, const real y, const real z) { return Vector2(x, y); }
		template<> Vector3 V<3>(const real x, const real y, const real z) { return Vector3(x, y, z); }

		template<> Vector<real, 1> V<1>(const Vector2 v2) { return Vector1(v2[0]); }
		template<> Vector<real, 2> V<2>(const Vector2 v2) { return v2; }
		template<> Vector<real, 3> V<3>(const Vector2 v2) { return Vector3(v2[0], v2[1], (real)0); }
		template<> Vector<real, 1> V<1>(const Vector3 v3) { return Vector1(v3[0]); }
		template<> Vector<real, 2> V<2>(const Vector3 v3) { return Vector2(v3[0], v3[1]); }
		template<> Vector<real, 3> V<3>(const Vector3 v3) { return v3; }

		template<> Vector1i Vi<1>(const int x, const int y, const int z, const int w) { return Vector1i(x); }
		template<> Vector2i Vi<2>(const int x, const int y, const int z, const int w) { return Vector2i(x, y); }
		template<> Vector3i Vi<3>(const int x, const int y, const int z, const int w) { return Vector3i(x, y, z); }
		template<> Vector4i Vi<4>(const int x, const int y, const int z, const int w) { return Vector4i(x, y, z, w); }
	}

}