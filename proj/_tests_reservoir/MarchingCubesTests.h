#pragma once

#include "Grid.h"
#include "Field.h"
#include "Mesh.h"
#include "MarchingCubes.h"
#include "Timer.h"
namespace Meso {


	template<class T, int d>
	void Test_Marching_Cubes_CPU(int times) {
		Typedef_VectorD(d);Timer timer;

		Vector<int, d> counts = Vector3i(100, 100, 100);
		VectorD domain_min = MathFunc::V<d>(0., 0., 0.);
		Grid<d> grid(counts, 0.01, domain_min, COLLOC);
		Field<T, d> field(grid, 0);

		Array<T> my_data(grid.DoF());
		VectorD center = MathFunc::V<d>(0.6, 0.6, 0.);
		grid.Exec_Nodes(
			[&](const VectorDi node) {
				VectorD pos = grid.Position(node);
				int index = grid.Index(node);
				field.Data()[index] = (pos - center).norm() - 0.6;
			}
		);
		auto m = std::make_shared<TriangleMesh<d>>();
		
		timer.Reset();
		for (size_t i = 0; i < times; i++)
		{
			Marching_Cubes<T, d>(field, m);
			double time = timer.Lap_Time(PhysicalUnits::s);
			double total_time = timer.Total_Time(PhysicalUnits::s);
			// Info("Used time: {:.2f}s/{:.2f}s, ETA {:.2f}s", time, total_time, total_time / (i + 1));
		}

		OBJFunc::Write_Mesh("./marching_cubes.obj", m);
		Pass("Test_Marching_Cubes Passed!");

		return;
	}


	template<class T, int d>
	void Test_Marching_Cubes_GPU(int times) {
		Typedef_VectorD(d);Timer timer;

		Vector<int, d> counts = Vector3i(100, 100, 100);
		VectorD domain_min = MathFunc::V<d>(0., 0., 0.);
		Grid<d> grid(counts, 0.01, domain_min, COLLOC);
		Field<T, d> field_host(grid, 0);
		Field<T, d, DEVICE> field(grid, 0);
		VectorD center = MathFunc::V<d>(0.6, 0.6, 0.);
		grid.Exec_Nodes(
			[&](const VectorDi node) {
				VectorD pos = grid.Position(node);
				int index = grid.Index(node);
				field_host.Data()[index] = (pos - center).norm() - 0.6;
			}
		);

		field = field_host;
		auto m = std::make_shared<TriangleMesh<d>>();

		timer.Reset();
		for (size_t i = 0; i < times; i++)
		{
			Marching_Cubes_GPU<T, d>(field, m);
			double time = timer.Lap_Time(PhysicalUnits::s);
			double total_time = timer.Total_Time(PhysicalUnits::s);
			// Info("Used time: {:.2f}s/{:.2f}s, ETA {:.2f}s", time, total_time, total_time / (i+1));
		}
		OBJFunc::Write_Mesh("./marching_cubes_GPU.obj", m);

		Pass("Test_Marching_Cubes[GPU] Passed!");

		return;
	}

	template<class T, int d>
	void Test_Marching_Cubes() {
		int times = 1;
		Test_Marching_Cubes_CPU<T, d>(times);
		Test_Marching_Cubes_GPU<T, d>(times);


	}
}