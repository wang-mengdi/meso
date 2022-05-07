#pragma once

#include "Grid.h"
#include "Field.h"
#include "Mesh.h"
#include "MarchingCubes.h"
namespace Meso {
	template<class T, int d>
	void Test_Marching_Cubes() {
		Typedef_VectorD(d);

		Vector<int, d> counts = Vector3i(12, 12, 12);
		VectorD domain_min = VectorFunc::V<d>(0., 0., 0.);
		Grid<d> grid(counts, 0.1, domain_min, COLLOC);

		Array<T> my_data(grid.DoF());
		VectorD center = VectorFunc::V<d>(0.6, 0.6, 0.);
		grid.Exec_Nodes(
			[&](const VectorDi node) {
				VectorD pos = grid.Position(node);
				int index = grid.Index(node);
				my_data[index] = (pos - center).norm() - 0.6;
			}
		);
		Field<T, d> field(grid, std::make_shared<Array<T>>(my_data));
		auto m = std::make_shared<TriangleMesh<d>>();
		Marching_Cubes<T, d>(field, m);
		OBJFunc::Write_Mesh("./marching_cubes.obj", m);
		return;
	}
}