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
		Array<VectorEi> elements;
		ArrayPtr<VectorD> normals; // normals in each faces

		// constructors(defalut: deep copy)
		SimplicialMesh(const ArrayPtr<VectorD> _vertices = nullptr) {vertices = _vertices == nullptr ? std::make_shared<Array<VectorD>>() : _vertices;}
		SimplicialMesh(const SimplicialMesh<d, e_d, side>& copy) { *this = copy; }

		SimplicialMesh<d, e_d, side>& operator=(const SimplicialMesh<d, e_d, side>& copy) {
			if (vertices == nullptr) { vertices = std::make_shared<Array<VectorD> >(); }
			*vertices = *(copy.vertices); 
			elements = copy.elements; 
			return *this; 
		}

		// attribute access
		static constexpr int Dim() { return d; }
		static constexpr int Face_Dim() { return e_d; }

		virtual Array<VectorD>& Vertices() { return *vertices; }
		virtual const Array<VectorD>& Vertices() const{ return *vertices; }

		virtual Array<VectorEi>& Elements() { return elements; }
		virtual const Array<VectorEi>& Elements() const{ return elements; }

		virtual const Array<VectorD>& Normals() const{ return *normals; }

		virtual void Clear() {
			if (vertices)vertices->clear(); 
			elements.clear();
			if (normals)normals->clear();
		}
	};

	template<int d, DataHolder side = DataHolder::HOST> class SegmentMesh : public SimplicialMesh<d, 2, side> {
		using Base = SimplicialMesh<d, 2,side>; Typedef_VectorD(d);
	public:
		SegmentMesh(const ArrayPtr<VectorD> _vertices = nullptr) : Base(_vertices) {}
	};

	template<int d, DataHolder side = DataHolder::HOST> class TriangleMesh : public SimplicialMesh<d, 3, side> {
		using Base = SimplicialMesh<d, 3, side>; Typedef_VectorD(d);
	public:
		TriangleMesh(const ArrayPtr<VectorD> _vertices = nullptr) : Base(_vertices) {}
	};

	template<int d> using SurfaceMesh = typename If<d == 2, SegmentMesh<2>, TriangleMesh<3> >::Type;
}