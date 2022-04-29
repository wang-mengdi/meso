#pragma once
namespace ThinShellAuxFunc {
	//p0-p1 is the shared edge, range: [-pi,+pi]
	inline real Dihedral_Angle(const Vector3& n0, const Vector3& n1, const Vector3& p0, const Vector3& p1) {
		real cosine = n0.dot(n1);
		Vector3 e = (p0 - p1).normalized();
		real sine = e.dot(n0.cross(n1));
		return atan2(sine, cosine);
	}

	inline int Third_Vertex(int p0, int p1, const Vector3i& vtxs) {
		for (int i = 0; i < 3; i++) {
			if (p0 != vtxs[i] && p1 != vtxs[i]) { return vtxs[i]; }
		}
		return -1; //shouldn't reach here
	}

	inline Vector2 Barycentric_Weights(const Vector3& p, const Vector3& p0, const Vector3& p1) {
		Vector3 e = p1 - p0;
		double t = e.dot(p - p0) / e.dot(e);
		return Vector2(1 - t, t);
	}

	//distance from one vertex to base, p0-p1 is the base, always positive
	inline real Distance(const Vector3& p, const  Vector3& p0, const  Vector3& p1) {
		return ((p - p0).cross(p - p1)).norm() / (p0 - p1).norm();
	}
};
