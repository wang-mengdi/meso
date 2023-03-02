//////////////////////////////////////////////////////////////////////////
// Test leveset on GPU and this is advanced test by senario examples
// Copyright (c) (2022-), Zhiqi Li, Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "LevelSet.h"
#include "Timer.h"
#include "GridEulerFunc.h"
#include "IOFunc.h"

namespace Meso {
	template<int d>
	void Fill_Eikonal_Error(Field<real, d>& error, const LevelSet<d>& levelset) {
		Typedef_VectorD(d);
		//calculate the residual of Eikonal equation by Rouy-Tourin
		error.Init(levelset.phi.grid);
		error.Calc_Nodes(
			[&](const VectorDi cell)->real {
				real phi0 = levelset.phi(cell);
				real sum = 0;
				for (int axis = 0; axis < d; axis++) {
					real diff = 0;
					for (int side = 0; side < 2; side++) {
						VectorDi nb_cell = Grid<d>::Neighbor_Node(cell, axis, side);
						if (!levelset.phi.grid.Valid(nb_cell)) continue;
						real phi1 = levelset.phi(nb_cell);
						//don't check interface cells since they may be wrong
						if (MathFunc::Sign(phi0) != MathFunc::Sign(phi1)) return 0;
						diff = std::max(diff, (phi0 - phi1) * MathFunc::Sign(phi0));
					}
					sum += diff * diff;
				}
				sum = sqrt(sum / (levelset.phi.grid.dx * levelset.phi.grid.dx));
				return sum - 1.0;
			}
		);
	}

	template<int d, DataHolder side>
	void Test_Fast_Marching(Grid<d> grid, const ImplicitManifold<d>& shape) {
		Typedef_VectorD(d);
		LevelSet<d> levelset(grid);
		//fill the levelset
		levelset.phi.Calc_Nodes(
			[&](const VectorDi cell)->real {
				VectorD pos = grid.Position(cell);
				real phi = shape.Phi(pos);
				if (std::fabs(phi) < 3 * grid.dx) return phi;
				else return phi * 2;
				//else return std::numeric_limits<real>::max();
				//else return -1;
			}
		);

		//fast marching
		Timer timer;
		levelset.Fast_Marching(-1);
		real fmm_time = timer.Lap_Time(PhysicalUnits::s);


		Field<real, d> fmm_error;
		Fill_Eikonal_Error(fmm_error, levelset);
		real max_eikonal_error = GridEulerFunc::Linf_Norm(fmm_error);

		LevelSet<d> analytical_levelset;
		analytical_levelset.Init(grid, shape);
		Field<real, d> distance_error(grid);
		distance_error = levelset.phi;
		distance_error -= analytical_levelset.phi;

		real max_distance_error = GridEulerFunc::Linf_Norm(distance_error);

		vtkNew<vtkStructuredGrid> vtk_grid;
		VTKFunc::VTS_Init_Grid(*vtk_grid, grid.Corner_Grid());
		VTKFunc::VTS_Add_Field(*vtk_grid, distance_error, "distance_error");
		VTKFunc::VTS_Add_Field(*vtk_grid, fmm_error, "fmm_error");

		std::stringstream filename;
		filename << "./fast_marching_"<<grid.Counts().transpose();
		VTKFunc::Write_VTS(*vtk_grid, filename.str() + ".vts");

		Pass("Fast Marching test passed for counts={} in {}s with eikonal linf error={}, distance linf error={}", grid.Counts(), fmm_time, max_eikonal_error, max_distance_error);
	}

	//Mengdi
	//See: Some Improvements of the Fast Marching Method
	template<int d>
	void Test_Fast_Marching(int scale) {
		Typedef_VectorD(d);
		real dx = 4.0 / scale;
		VectorD domain_min = MathFunc::V<d>(-2, -2, -2);
		VectorDi counts = MathFunc::Vi<d>(scale, scale, scale);
		Grid<d> grid(counts, dx, domain_min, CENTER);
		Sphere<d> sphere(VectorD::Zero(), 0.5);
		Test_Fast_Marching<d, HOST>(grid, sphere);
	}
}