#pragma once

#include "Grid.h"
#include "Field.h"
#include "Mesh.h"
#include "MarchingCubes.h"
namespace Meso {
	template<class T, int d>
	void Test_Marching_Cubes() {
		Typedef_VectorD(d);

		Vector<int, d> counts = VectorDi(100, 100, 100);
		VectorD domain_min = VectorD(0., 0., 0.);
		Grid<d> grid(counts, 0.01, domain_min, COLLOC);
		Field<T, d> field(grid);

		VectorD center = VectorD(0.6, 0.6, 0.);
		grid.Exec_Nodes(
			[&](const VectorDi node) {
				VectorD pos = grid.Position(node);
				int index = grid.Index(node);
				field.Data()[index] = (pos - center).norm() - 0.6;
			}
		);
		auto m = std::make_shared<TriangleMesh<d>>();
		Marching_Cubes<T, d>(field, m);
		OBJFunc::Write_Mesh("./marching_cubes.obj", m);

		Pass("Test_Marching_Cubes Passed!");

		return;
	}


	template<class T, int d>
	void Test_Marching_Cubes_GPU() {
		Typedef_VectorD(d);

		Vector<int, d> counts = VectorDi(100, 100, 100);
		VectorD domain_min = VectorD(0., 0., 0.);
		Grid<d> grid(counts, 0.01, domain_min, COLLOC);
		Field<T, d, DEVICE> field(grid);

		VectorD center = VectorD(0.6, 0.6, 0.);
		grid.Exec_Nodes(
			[&](const VectorDi node) {
				VectorD pos = grid.Position(node);
				int index = grid.Index(node);
				field.Data()[index] = (pos - center).norm() - 0.6;
			}
		);
		auto m = std::make_shared<TriangleMesh<d>>();
		Marching_Cubes_GPU<T, d>(field, m);
		OBJFunc::Write_Mesh("./marching_cubes_GPU.obj", m);

		Pass("Test_Marching_Cubes[GPU] Passed!");

		return;
	}
}