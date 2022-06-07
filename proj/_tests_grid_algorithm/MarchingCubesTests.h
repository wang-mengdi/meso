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

	template<class T, DataHolder side, int d>
	void Test_Marching_Cubes_Unit(const Field<T, d, side>& field, int times, bool verbose) {
		VertexMatrix<T, d> vertices; ElementMatrix<d> faces; Timer timer; timer.Reset();
		for (size_t i = 0; i < times; i++)
		{
			Marching_Cubes<T, d, side>(vertices, faces, field);
			double time = timer.Lap_Time(PhysicalUnits::s);
			double total_time = timer.Total_Time(PhysicalUnits::s);
			if (verbose) Info("Used time: {:.2f}s/{:.2f}s, ETA {:.2f}s", time, total_time, total_time / (i + 1));
		}
		std::stringstream filename;
		filename << "./marching_" << ((d == 2) ? "square" : "cubes") << "_" << ((side == HOST) ? "CPU" : "GPU") << ".obj";
		//sprintf(filename, "./marching_%s_%s.obj", (d == 2) ? "square" : "cubes", (side == HOST) ? "CPU" : "GPU");
		OBJFunc::Write_OBJ(filename.str(), vertices, faces);
		Pass("Test_Marching_Cubes[{}] Passed with {} vertices and {} elements", (side == HOST) ? "CPU" : "GPU", vertices.rows(), faces.rows());
	}

	template<class T>
	void Test_Marching_Cubes(int times, bool verbose) {

		// Marching Cubes
		{
			Grid<3> grid(Vector3i(100, 100, 100), 0.01, MathFunc::V<3>(0., 0., 0.), COLLOC);
			Field<T, 3> field_host(grid, 0); Vector3 center = MathFunc::V<3>(0.6, 0.6, 0.);

			grid.Exec_Nodes(
				[&](const Vector3i node) {
					Vector3 pos = grid.Position(node);
					field_host(node) = (pos - center).norm() - 0.6;
				}
			);

			Test_Marching_Cubes_Unit<T, DEVICE, 3>(Field<T, 3, DEVICE>(field_host), times, verbose);
			Test_Marching_Cubes_Unit<T, HOST, 3>(field_host, times, verbose);

		}

		// Marching Square
		{
			Grid<2> grid(Vector2i(100, 100), 0.01, MathFunc::V<2>(0., 0.), COLLOC);
			Field<T, 2> field_host(grid, 0); Vector2 center = grid.Center();
			grid.Exec_Nodes(
				[&](const Vector2i node) {
					Vector2 pos = (grid.Position(node) - center);
					field_host(node) = 1000. * pos[1] * pos[1] - 300. * pos[0] * pos[0] * sin(100. * pos[0] * pos[0]);
				}
			);
			Test_Marching_Cubes_Unit<T, HOST, 2>(field_host, times, verbose);

		}
	}

}