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

		Array<std::string> Split_String(const std::string& s, const std::string& delimiters)
		{
			//https://stackoverflow.com/questions/26328793/how-to-split-string-with-delimiter-using-c
			Array<std::string> tokens; tokens.clear();
			std::string::size_type lastPos = s.find_first_not_of(delimiters, 0);
			std::string::size_type pos = s.find_first_of(delimiters, lastPos);
			while (std::string::npos != pos || std::string::npos != lastPos) {
				tokens.push_back(s.substr(lastPos, pos - lastPos));
				lastPos = s.find_first_not_of(delimiters, pos);
				pos = s.find_first_of(delimiters, lastPos);
			}
			return tokens;
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

		Vector1 Orthogonal_Vector(const Vector1& v) { return v; }
		Vector2 Orthogonal_Vector(const Vector2& v) { return { -v.y(),v.x() }; }
		Vector3 Orthogonal_Vector(const Vector3& v)
		{
			real abs_x = abs(v.x()), abs_y = abs(v.y()), abs_z = abs(v.z());
			if (abs_x < abs_y) return abs_x < abs_z ? Vector3((real)0, v.z(), -v.y()) : Vector3(v.y(), -v.x(), (real)0);
			else return abs_y < abs_z ? Vector3(-v.z(), (real)0, v.x()) : Vector3(v.y(), -v.x(), (real)0);
		}
		Vector4 Orthogonal_Vector(const Vector4& v)
		{
			Vector4 n = Vector4::Zero(); int min_axis = 0; real min_abs = abs(v[0]); for (int i = 1; i < 4; i++)if (abs(v[i]) < min_abs) { min_abs = abs(v[i]); min_axis = i; }
			Vector3 v3; {int c = 0; for (int i = 0; i < 4; i++) { if (i == min_axis)continue; v3[c++] = v[i]; }}
			v3 = Orthogonal_Vector(v3); {int c = 0; for (int i = 0; i < 4; i++) { if (i == min_axis)continue; n[i] = v3[c++]; }}return n;
		}

		Vector1 Cross(const Vector2& v1, const Vector2& v2) { Vector1 v; v[0] = (v1[0] * v2[1] - v1[1] * v2[0]); return v; }
		Vector2 Cross(const Vector1& v1, const Vector2& v2) { return Orthogonal_Vector(v2) * v1[0]; }
		Vector2 Cross(const real& v1, const Vector2& v2) { return Orthogonal_Vector(v2) * v1; }

		real Angle_Between(const Vector2& v1, const Vector2& v2)
		{
			real c = v1.dot(v2); return acos(c);
		}

		real Angle_Between(const Vector3& v1, const Vector3& v2)
		{
			real c = v1.dot(v2); return acos(c);
		}

		real Angle_From_To(const Vector2& v1, const Vector2& v2)
		{
			real s = Cross(v1,v2)[0];
			real c = v1.dot(v2);
			return atan2(s, c);
		}

		real Angle_From_To(const Vector3& v1, const Vector3& v2)
		{
			real s = (v1.cross(v2)).norm();
			real c = v1.dot(v2);
			return atan2(s, c);
		}
	}

}