//////////////////////////////////////////////////////////////////////////
// Mesh functions
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#include <numeric>
#include <iostream>
#include "ImplicitGeometry.h"
#include "Constants.h"
#include "Hashtable.h"
#include "Mesh.h"
#include "MeshFunc.h"

namespace Meso {
	namespace MeshFunc {
		////////////////////////////////////////////////////////////////////////////////////////////////////
		////Normals
		Vector2 Normal(const Vector2& v1, const Vector2& v2)
		{
			Vector2 v12 = (v2 - v1).normalized(); return Vector2(v12[1], -v12[0]);
		}

		Vector3 Normal(const Vector3& p1, const Vector3& p2, const Vector3& p3)
		{
			return (p2 - p1).cross(p3 - p1).normalized();
		}

		int Element_Edges(const Vector2i& vtx, Array<Vector2i>& edges)
		{
			edges.resize(1);
			edges[0] = vtx; return 1;
		}

		int Element_Edges(const Vector3i& vtx, Array<Vector2i>& edges)
		{
			edges.resize(3);
			edges[0] = Vector2i(vtx[0], vtx[1]); edges[1] = Vector2i(vtx[1], vtx[2]); edges[2] = Vector2i(vtx[2], vtx[0]); return 3;
		}
	};
}