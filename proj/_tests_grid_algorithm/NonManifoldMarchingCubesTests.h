//////////////////////////////////////////////////////////////////////////
// Test marching cubes algorithm
// Copyright (c) (2023-), Fan Feng
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "NonManifoldMarchingCubes.h"
#include "Timer.h"
#include "IOFunc.h"
#include "ImplicitManifold.h"

namespace Meso {

	template<class T, DataHolder side, int d>
	void Test_Non_Manifold_Marching_Cubes_Unit(const Field<CellType, d, side>& label, const Field<T, d, side>& field, int test_case, int times, bool verbose) {
		VertexMatrix<T, d> vertices; ElementMatrix<d> faces; Timer timer; timer.Reset();
		for (size_t i = 0; i < times; i++)
		{
			Non_Manifold_Marching_Cubes<T, d, side>(vertices, faces, label,field);
			double time = timer.Lap_Time(PhysicalUnits::s);
			double total_time = timer.Total_Time(PhysicalUnits::s);
			if (verbose) Info("Used time: {:.2f}s/{:.2f}s, ETA {:.2f}s", time, total_time, total_time / (i + 1));
		}
		std::stringstream filename;
		filename << "./non_manifold_marching_" << ((d == 2) ? "square" : "cubes") << "_" << ((side == HOST) ? "CPU" : "GPU") << test_case;
		OBJFunc::Write_OBJ(filename.str()+".obj", vertices, faces);
		vtkNew<vtkStructuredGrid> vtk_grid;
		VTKFunc::VTS_Init_Grid(*vtk_grid,field.grid);
		VTKFunc::VTS_Add_Field(*vtk_grid,field,"value");
		VTKFunc::VTS_Add_Field(*vtk_grid, label, "label");
		VTKFunc::Write_VTS(*vtk_grid, filename.str() +".vts");
		Pass("Test_Marching_{}[{}] Passed with {} vertices and {} elements", (d == 2) ? "Square" : "Cubes", (side == HOST) ? "CPU" : "GPU", vertices.rows(), faces.rows());
	}

	template<class T>
	void Test_Non_Manifold_Marching_Cubes1(int times, bool verbose) {

		// Marching Cubes
		/*{
			Grid<3> grid(Vector3i(100, 100, 100), 0.01, MathFunc::V<3>(0., 0., 0.), CORNER);
			Field<T, 3> field_host(grid, 0); Vector3 center = MathFunc::V<3>(0.6, 0.6, 0.);

			grid.Exec_Nodes(
				[&](const Vector3i node) {
					Vector3 pos = grid.Position(node);
					field_host(node) = (pos - center).norm() - 0.6;
				}
			);

			Test_Marching_Cubes_Unit<T, DEVICE, 3>(Field<T, 3, DEVICE>(field_host), times, verbose);
			Test_Marching_Cubes_Unit<T, HOST, 3>(field_host, times, verbose);
		}*/

		// Marching Square
		{
			Grid<2> grid(Vector2i(104, 104), 0.01, MathFunc::V<2>(0., 0.), CORNER);
			Field<T, 2> field_host(grid, 0); Vector2 center = grid.Center();
			Field<CellType, 2> label_host(grid, 0);
			grid.Exec_Nodes(
				[&](const Vector2i node) {
					Vector2 pos = (grid.Position(node) - center);
					real value= 1000. * pos[1] * pos[1] - 300. * pos[0] * pos[0] * sin(100. * pos[0] * pos[0]);
					field_host(node) = value;
					if (value > 0) { label_host(node) = 0; }
					else { field_host(node) = -value; label_host(node) = 1; }
				}
			);
			Test_Non_Manifold_Marching_Cubes_Unit<T, HOST, 2>(label_host, field_host, 1, times, verbose);
		}
	}

	template<class T>
	void Test_Non_Manifold_Marching_Cubes2(int times, bool verbose) {

		// Marching Cubes
		/*{
			Grid<3> grid(Vector3i(100, 100, 100), 0.01, MathFunc::V<3>(0., 0., 0.), CORNER);
			Field<T, 3> field_host(grid, 0); Vector3 center = MathFunc::V<3>(0.6, 0.6, 0.);

			grid.Exec_Nodes(
				[&](const Vector3i node) {
					Vector3 pos = grid.Position(node);
					field_host(node) = (pos - center).norm() - 0.6;
				}
			);

			Test_Marching_Cubes_Unit<T, DEVICE, 3>(Field<T, 3, DEVICE>(field_host), times, verbose);
			Test_Marching_Cubes_Unit<T, HOST, 3>(field_host, times, verbose);
		}*/

		// Marching Square
		{
			Grid<2> grid(Vector2i(104, 104), 0.01, MathFunc::V<2>(0., 0.), CORNER);
			Field<T, 2> field_host(grid, 0); Vector2 center1 = Vector2(0.5,0.3); Vector2 center2 = Vector2(0.33,0.67); Vector2 center3 = Vector2(0.67,0.67);
			Field<CellType, 2> label_host(grid, 0);
			PlaneShape plane12 = PlaneShape<2>(0.5*(center1+center2),center1-center2);
			PlaneShape plane13 = PlaneShape<2>(0.5 * (center1 + center3), center1 - center3);
			PlaneShape plane23 = PlaneShape<2>(0.5 * (center3 + center2), center3 - center2);

			grid.Exec_Nodes(
				[&](const Vector2i node) {
					Vector2 pos = grid.Position(node);
					real dist1 = -(pos - center1).norm();
					real dist2= -(pos - center2).norm();
					real dist3= -(pos - center3).norm();
					real max_value = std::max({ dist1, dist2, dist3 });
					if (max_value == dist1) {
						label_host(node) = 0;
						field_host(node) = -min(abs(plane12.Phi(pos)),abs(plane13.Phi(pos)));
					}
					else if (max_value == dist2) {
						label_host(node) = 1;
						field_host(node) = -min(abs(plane12.Phi(pos)), abs(plane23.Phi(pos)));
					}
					else {
						label_host(node) = 2;
						field_host(node) = -min(abs(plane13.Phi(pos)), abs(plane23.Phi(pos)));
					}
				}
			);
			Test_Non_Manifold_Marching_Cubes_Unit<T, HOST, 2>(label_host, field_host, 2, times, verbose);
		}
	}
}