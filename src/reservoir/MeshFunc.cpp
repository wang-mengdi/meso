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

		template<> Vector<real, 2> Normal<2>(const ArrayF<Vector<real, 2>, 2>& v) { return Normal(v[0], v[1]); }
		template<> Vector<real, 3> Normal<3>(const ArrayF<Vector<real, 3>, 3>& v) { return Normal(v[0], v[1], v[2]); }

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

		//////////////////////////////////////////////////////////////////////////
		////Mesh edges
		template<int d, int e_d> void Get_Edges(const SimplicialMesh<d, e_d>& mesh, Array<Vector2i>& edges)
		{
			Typedef_VectorEi(e_d);
			Hashset<Vector2i> edge_hashset; Array<Vector2i> element_edges;
			for (const VectorEi& vtx : mesh.elements) {
				int n = Element_Edges(vtx, element_edges);
				for (int i = 0; i < n; i++) {
					edge_hashset.insert(Sorted(element_edges[i]));
				}
			}
			for (const auto& edge : edge_hashset) { edges.push_back(edge); }
		}

#define Inst_Helper(d,e_d) \
	template void Get_Edges(const SimplicialMesh<d,e_d>&,Array<Vector2i>&);
		Inst_Helper(2, 2); Inst_Helper(2, 3); Inst_Helper(3, 2); Inst_Helper(3, 3);
#undef Inst_Helper

		//////////////////////////////////////////////////////////////////////////
		////Mesh initialization
		template<class T_MESH> void Initialize_Herring_Bone_Mesh(const int m, const int n, const real dx, T_MESH* mesh, int axis_0, int axis_1)
		{
			using VectorD = Vector<real, T_MESH::Dim()>;
			mesh->elements.resize(2 * (m - 1) * (n - 1)); int t = 0;
			for (int i = 1; i <= m - 1; i++)for (int j = 1; j <= n - 1; j++) { // counterclockwise node ordering
				if (i % 2) { mesh->elements[t++] = Vector3i(i + m * (j - 1), i + 1 + m * (j - 1), i + m * j); mesh->elements[t++] = Vector3i(i + 1 + m * (j - 1), i + 1 + m * j, i + m * j); }
				else { mesh->elements[t++] = Vector3i(i + m * (j - 1), i + 1 + m * (j - 1), i + 1 + m * j); mesh->elements[t++] = Vector3i(i + m * (j - 1), i + 1 + m * j, i + m * j); }
			}
			for (int i = 0; i < mesh->elements.size(); i++) { mesh->elements[i] -= Vector3i::Ones(); }
			for (int j = 0; j < n; j++)for (int i = 0; i < m; i++) { VectorD pos = VectorD::Zero(); pos[axis_0] = (real)i * dx; pos[axis_1] = (real)j * dx; mesh->Vertices().push_back(pos); }
		}

#define Inst_Helper(T1) \
	template void Initialize_Herring_Bone_Mesh<T1>(const int,const int,const real,T1*,int,int);
		Inst_Helper(TriangleMesh<2>);
		Inst_Helper(TriangleMesh<3>);
#undef Inst_Helper
	};
}