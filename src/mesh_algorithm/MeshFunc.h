//////////////////////////////////////////////////////////////////////////
// Mesh functions
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "Common.h"
#include "Mesh.h"

namespace Meso {
    namespace MeshFunc {
        ////Normals
        Vector2 Normal(const Vector2& v1, const Vector2& v2);
        Vector3 Normal(const Vector3& v1, const Vector3& v2, const Vector3& v3);
        template<int d> Vector<real, d> Normal(const ArrayF<Vector<real, d>, d>& v) {
            if constexpr (d == 2) {
                return Normal(v[0], v[1]);
            }
            else if constexpr (d == 3) {
                return Normal(v[0], v[1], v[2]);
            }
            else {
                Error("Dimension not supported.");
            }
        }

        ////Element edges
        int Element_Edges(const Vector2i& v, Array<Vector2i>& edges);
        int Element_Edges(const Vector3i& v, Array<Vector2i>& edges);

        ////Mesh edges
        template<int d, int e_d> void Get_Edges(const SimplicialMesh<d, e_d>& mesh, Array<Vector2i>& edges) {
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

        ////Mesh initialization with analytical shapes
        template<class MeshType> void Initialize_Herring_Bone_Mesh(const int m, const int n, const real dx, MeshType* mesh, int axis_0 = 0, int axis_1 = 1) {
            using VectorD = Vector<real, MeshType::Dim()>;
            mesh->elements.resize(2 * (m - 1) * (n - 1)); int t = 0;
            for (int i = 1; i <= m - 1; i++)for (int j = 1; j <= n - 1; j++) { // counterclockwise node ordering
                if (i % 2) { mesh->elements[t++] = Vector3i(i + m * (j - 1), i + 1 + m * (j - 1), i + m * j); mesh->elements[t++] = Vector3i(i + 1 + m * (j - 1), i + 1 + m * j, i + m * j); }
                else { mesh->elements[t++] = Vector3i(i + m * (j - 1), i + 1 + m * (j - 1), i + 1 + m * j); mesh->elements[t++] = Vector3i(i + m * (j - 1), i + 1 + m * j, i + m * j); }
            }
            for (int i = 0; i < mesh->elements.size(); i++) { mesh->elements[i] -= Vector3i::Ones(); }
            for (int j = 0; j < n; j++)for (int i = 0; i < m; i++) { VectorD pos = VectorD::Zero(); pos[axis_0] = (real)i * dx; pos[axis_1] = (real)j * dx; mesh->Vertices().push_back(pos); }
        }
    };
}