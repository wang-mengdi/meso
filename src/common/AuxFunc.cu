#include "AuxFunc.h"

namespace Meso {
	
	namespace IOFunc {
		void Create_Directory(const bf::path path)
		{
			//recursively
			try {
				if (!boost::filesystem::exists(path))
					boost::filesystem::create_directories(path);
			}
			catch (std::exception& e) { // Not using fs::filesystem_error since std::bad_alloc can throw too.
				std::cout << e.what() << std::endl;
			}
		}
		std::string To_String_Simple(const bool& a) {
			return std::to_string((int)a);
		}
	}

	namespace VectorFunc {
		template<> __host__ __device__ Vector1 V<1>(const real x, const real y, const real z) { return Vector1(x); }
		template<> __host__ __device__ Vector2 V<2>(const real x, const real y, const real z) { return Vector2(x, y); }
		template<> __host__ __device__ Vector3 V<3>(const real x, const real y, const real z) { return Vector3(x, y, z); }

		template<> __host__ __device__ Vector<real, 1> V<1>(const Vector2 v2) { return Vector1(v2[0]); }
		template<> __host__ __device__ Vector<real, 2> V<2>(const Vector2 v2) { return v2; }
		template<> __host__ __device__ Vector<real, 3> V<3>(const Vector2 v2) { return Vector3(v2[0], v2[1], (real)0); }
		template<> __host__ __device__ Vector<real, 1> V<1>(const Vector3 v3) { return Vector1(v3[0]); }
		template<> __host__ __device__ Vector<real, 2> V<2>(const Vector3 v3) { return Vector2(v3[0], v3[1]); }
		template<> __host__ __device__ Vector<real, 3> V<3>(const Vector3 v3) { return v3; }
	}

}