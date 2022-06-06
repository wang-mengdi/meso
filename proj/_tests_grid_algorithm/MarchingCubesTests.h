//////////////////////////////////////////////////////////////////////////
// Test marching cubes algorithm
// Copyright (c) (2022-), Yunquan Gu
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once

#include "MarchingCubes.h"
#include "Timer.h"
#include "IOFunc.h"

namespace Meso {


	template<class T, int d, DataHolder side>
	void Test_Marching_Cubes(int times, bool verbose) {
		Typedef_VectorD(d);

		Timer timer;
		Vector<int, d> counts = Vector3i(100, 100, 100);
		VectorD domain_min = MathFunc::V<d>(0., 0., 0.);
		Grid<d> grid(counts, 0.01, domain_min, COLLOC);
		Field<T, d> field_host(grid, 0);

		Array<T> my_data(grid.DoF());
		VectorD center = MathFunc::V<d>(0.6, 0.6, 0.);
		grid.Exec_Nodes(
			[&](const VectorDi node) {
				VectorD pos = grid.Position(node);
				int index = grid.Index(node);
				field_host.Data()[index] = (pos - center).norm() - 0.6;
			}
		);

		Field<T, d, side> field = field_host;

		VertexMatrix<T, d> vertices;
		ElementMatrix<3> faces;
		
		timer.Reset();
		for (size_t i = 0; i < times; i++)
		{
			Marching_Cubes<T, d, side>(vertices, faces, field);
			double time = timer.Lap_Time(PhysicalUnits::s);
			double total_time = timer.Total_Time(PhysicalUnits::s);
			if (verbose) Info("Used time: {:.2f}s/{:.2f}s, ETA {:.2f}s", time, total_time, total_time / (i + 1));
		}

		//TriangleMesh<d> mesh;
		//mesh.vertices.resize(vertices.rows());
		//mesh.elements.resize(faces.rows());
		//for (int i = 0; i < vertices.rows(); i++) {
		//	for (int j = 0; j < d; j++) {
		//		mesh.vertices[i][j] = vertices(i, j);
		//	}
		//}
		//for (int i = 0; i < faces.rows(); i++) {
		//	for (int j = 0; j < 3; j++) {
		//		mesh.elements[i][j] = faces(i, j);
		//	}
		//}

		std::string output_name;
		if (side == DEVICE) output_name = "./marching_cubes_GPU.obj";
		else output_name = "./marching_cubes_CPU.obj";
		OBJFunc::Write_OBJ(output_name, vertices, faces);

		if (side == HOST) Pass("Test_Marching_Cubes[CPU] Passed with {} vertices and {} elements", vertices.rows(), faces.rows());
		else Pass("Test_Marching_Cubes[GPU] Passed with {} vertices and {} elements", vertices.rows(), faces.rows());
	}

	//template<class T, int d>
	//void Test_Marching_Cubes() {
	//	bool verbose = false; int times = 1;
	//	Test_Marching_Cubes_CPU<T, d>(times, verbose);
	//	Test_Marching_Cubes_GPU<T, d>(times, verbose);
	//}
}