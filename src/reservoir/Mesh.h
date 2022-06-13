//////////////////////////////////////////////////////////////////////////
// Mesh
// Copyright (c) (2022-) Bo zhu, Yunquan Gu, Mengdi Wang, Fan Feng
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include <Common.h>

namespace Meso {
	template<class T, int d> using VertexMatrix = Eigen::Matrix<T, Eigen::Dynamic, d, Eigen::RowMajor>; //each row is the coordinate of a triangle
	template<int e_d> using ElementMatrix = Eigen::Matrix<int, Eigen::Dynamic, e_d, Eigen::RowMajor>; //each row contains the indices of one element
	using TriangleElementMatrix = Eigen::Matrix<int, Eigen::Dynamic, 3, Eigen::RowMajor>;

	template<class T, int d, int e_d, DataHolder side = DataHolder::HOST>
	class SimplicialMesh {
		Typedef_VectorD(d);
		Typedef_VectorEi(e_d);
	public:
		VertexMatrix<T, d> vertices;
		ElementMatrix<e_d> elements;

		// constructors(defalut: deep copy)
		SimplicialMesh() {}
		SimplicialMesh(VertexMatrix<T, d> _vertices, ElementMatrix<e_d> _elements) { vertices = _vertices; elements = _elements; }
		SimplicialMesh(const SimplicialMesh<T, d, e_d, side>& copy) { *this = copy; }

		SimplicialMesh<T, d, e_d, side>& operator=(const SimplicialMesh<T, d, e_d, side>& copy) {
			vertices = copy.vertices;
			elements = copy.elements;
			return *this;
		}

		// attribute access
		static constexpr int Dim() { return d; }
		static constexpr int Face_Dim() { return e_d; }

		constexpr VertexMatrix<T, d>& Vertices() { return vertices; }
		constexpr const VertexMatrix<T, d>& Vertices() const { return vertices; }

		constexpr Vector<T, d>& Vertex(int i) { return vertices.row(i); }
		constexpr const Vector<T, d>& Vertex(int i) const { return vertices.row(i); }

		constexpr ElementMatrix<e_d>& Elements() { return elements; }
		constexpr const ElementMatrix<e_d>& Elements() const { return elements; }

		constexpr Vector<int, e_d>& Element(int i) { return elements.row(i); }
		constexpr const Vector<int,e_d>& Element(int i) const { return elements.row(i); }
	};

	template<int d, DataHolder side = DataHolder::HOST> using SegmentMesh = SimplicialMesh <real, d, 2, side>;
	template<int d, DataHolder side = DataHolder::HOST> using TriangleMesh = SimplicialMesh <real, d, 3, side>;

	template<int d> using SurfaceMesh = typename If<d == 2, SegmentMesh<2>, TriangleMesh<3> >::Type;
}