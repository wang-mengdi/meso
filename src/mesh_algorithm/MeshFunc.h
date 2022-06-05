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

        ////Element edges
        template<int d> int Element_Edges(const Vector<int, d>& vtx, Array<Vector2i>& edges) {
            if constexpr (d == 2) {
                edges.resize(1);
                edges[0] = vtx; return 1;
            }
            else if constexpr (d == 3) {
                edges.resize(3);
                edges[0] = Vector2i(vtx[0], vtx[1]); edges[1] = Vector2i(vtx[1], vtx[2]); edges[2] = Vector2i(vtx[2], vtx[0]); return 3;
            }
            else {
                Error("Dimension not supported");
            }
        }

        ////Mesh edges
        template<int d, int e_d> void Get_Edges(const SimplicialMesh<d, e_d>& mesh, Array<Vector2i>& edges) {
            Typedef_VectorEi(e_d);
            Hashset<Vector2i> edge_hashset; Array<Vector2i> element_edges;
            for (const VectorEi& vtx : mesh.elements) {
                int n = Element_Edges<e_d>(vtx, element_edges);
                for (int i = 0; i < n; i++) {
                    edge_hashset.insert(MathFunc::Sorted<2>(element_edges[i]));
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