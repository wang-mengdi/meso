//////////////////////////////////////////////////////////////////////////
// Mesh
// Copyright (c) (2022-) Bo zhu, Yunquan Gu
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once

#include <Common.h>

namespace Meso {
	template<int d, int e_d, DataHolder side = DataHolder::HOST>
	class SimplicialMesh {
		Typedef_VectorD(d);Typedef_VectorEi(e_d);
	public:
		ArrayPtr<VectorD> vertices;
		Array<VectorEi> faces;
		ArrayPtr<VectorD> normals; // normals in each faces

		// constructors(defalut: deep copy)
		SimplicialMesh() {}
		//SimplicialMesh(const ArrayPtr<VectorD> _vertices = nullptr) {
			//vertices = _vertices == nullptr ? std::make_shared<Array<VectorD>>() : _vertices;
		//}
		//SimplicialMesh(const SimplicialMesh<d, e_i, slide>& copy) { *this = copy; }
		//SimplicialMesh<d, e_i, slide>& opertor = (const SimplicialMesh<d, slide>&copy){ *this = copy; }

		// attribute access
		static constexpr int Dim() { return d; }
		static constexpr int Face_Dim() { return e_d; }

		virtual const Array<VectorD>& Vertices() { return *vertices; }
		virtual const Array<VectorEi>& Faces() { return faces; }
		virtual const Array<VectorD>& Normals() { return *normals; }

		virtual void Clear() {
			if (vertices)vertices->clear(); faces.clear();
			if (normals)normals->clear();
		}

	};

	template<int d> class TriangleMesh : public SimplicialMesh<d, 3> {
		using Base = SimplicialMesh<d, 3>; Typedef_VectorD(d); Typedef_VectorEi(3);
	public:
		using Base::vertices; using Base::faces; using Base::normals;
		TriangleMesh() { this->vertices = std::make_shared<Array<Vector3>>(); }
		/*TriangleMesh(const ArrayPtr<VectorD> _vertices = nullptr) :Base(_vertices) {}*/
	};
}