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
        template<int d> Vector<real, d> Normal(const ArrayF<Vector<real, d>, d>& v);

        ////Element edges
        int Element_Edges(const Vector2i& v, Array<Vector2i>& edges);
        int Element_Edges(const Vector3i& v, Array<Vector2i>& edges);

        ////Mesh edges
        template<int d, int e_d> void Get_Edges(const SimplicialMesh<d, e_d>& mesh, Array<Vector2i>& edges);

        ////Mesh initialization with analytical shapes
        template<class T_MESH> void Initialize_Herring_Bone_Mesh(const int m, const int n, const real dx, T_MESH* mesh, int axis_0 = 0, int axis_1 = 1);
    };
}