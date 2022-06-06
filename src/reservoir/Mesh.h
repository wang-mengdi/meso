//////////////////////////////////////////////////////////////////////////
// Mesh
// Copyright (c) (2022-) Bo zhu, Yunquan Gu
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include <Common.h>

namespace Meso {
	template<class T, int d> using VertexMatrix = Eigen::Matrix<T, Eigen::Dynamic, d, Eigen::RowMajor>;
	template<int d> using ElementMatrix = Eigen::Matrix<int, Eigen::Dynamic, d, Eigen::RowMajor>;

	template<int d, int e_d, DataHolder side = DataHolder::HOST>
	class SimplicialMesh {
		Typedef_VectorD(d);
		Typedef_VectorEi(e_d);
	public:
		Array<VectorD, side> vertices;
		Array<VectorEi, side> elements;
		//Array<VectorD, side> normals; // normals in each faces

		// constructors(defalut: deep copy)
		SimplicialMesh() {}
		SimplicialMesh(const ArrayPtr<VectorD>& _vertices) { vertices = _vertices; }
		SimplicialMesh(const SimplicialMesh<d, e_d, side>& copy) { *this = copy; }

		SimplicialMesh<d, e_d, side>& operator=(const SimplicialMesh<d, e_d, side>& copy) {
			vertices = (copy.vertices);
			elements = copy.elements;
			//normals = copy.normals;
			return *this;
		}

		// attribute access
		static constexpr int Dim() { return d; }
		static constexpr int Face_Dim() { return e_d; }

		constexpr Array<VectorD, side>& Vertices() { return vertices; }
		constexpr const Array<VectorD, side>& Vertices() const { return vertices; }

		constexpr Array<VectorEi, side>& Elements() { return elements; }
		constexpr const Array<VectorEi, side>& Elements() const { return elements; }

		//constexpr const Array<VectorD, side>& Normals() const { return normals; }

		void Clear() {
			vertices.clear(); 
			elements.clear();
			//normals.clear();
		}
	};

	template<int d, DataHolder side = DataHolder::HOST> using SegmentMesh = SimplicialMesh<d, 2, side>;
	template<int d, DataHolder side = DataHolder::HOST> using TriangleMesh = SimplicialMesh<d, 3, side>;

	template<int d> using SurfaceMesh = typename If<d == 2, SegmentMesh<2>, TriangleMesh<3> >::Type;
}