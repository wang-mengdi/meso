//////////////////////////////////////////////////////////////////////////
// Test leveset on GPU and this is advanced test by senario examples
// Copyright (c) (2022-), Zhiqi Li, Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "LevelSet.h"
#include "Timer.h"
#include "GridEulerFunc.h"

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
	void Test_Fast_Marching(Grid<d> grid, const ImplicitGeometry<d>& shape) {
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
		levelset.phi -= analytical_levelset.phi;

		real max_distance_error = GridEulerFunc::Linf_Norm(levelset.phi);

		//real max_err = fmm_error.Max_Abs();

		////for debugging
		//real max_err = -1;
		//VectorDi max_err_cell;
		//fmm_error.Iterate_Nodes(
		//	[&](const VectorDi cell) {
		//		if (std::fabs(fmm_error(cell)) > max_err) {
		//			max_err = std::fabs(fmm_error(cell));
		//			max_err_cell = cell;
		//		}
		//	}
		//);
		//Info("max error {} at cell {}", max_err, max_err_cell);

		Pass("Fast Marching test passed for counts={} in {}s with eikonal linf error={} and distance linf error={}", grid.Counts(), fmm_time, max_eikonal_error, max_distance_error);

		//real eps = 1;
		//real eps = sqrt(std::numeric_limits<real>::epsilon());
		//if (max_distance_error > eps) Error("Fast Marching for counts={} failed with max error={}", counts, max_err);
		//else Pass("Fast Marching passed for counts={} in {}s with max error={}", counts, fmm_time, max_err);
	}

	////Mengdi
	//template<int d, DataHolder side>
	//void Test_Fast_Marching(const Vector<int, d> counts) {
	//	Typedef_VectorD(d);
	//	VectorD domain_min = MathFunc::V<d>(-0.9, -1.2, 3);
	//	Grid<d> grid(counts, 0.01, domain_min, MAC);
	//	LevelSet<d> levelset(grid);
	//	VectorD domain_max = grid.Domain_Max();
	//	VectorD domain_len = domain_max - domain_min;
	//	real min_side = domain_len.minCoeff();
	//	//Sphere<d> sphere(domain_min + domain_len * 0.2, min_side * 0.3);
	//	Sphere<d> sphere(domain_min + domain_len * 0.5, min_side * 0.3);

	//	Test_Fast_Marching<d, side>(grid, sphere);
	//}


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

	///// Here, we test the fast marching method
	//template<int d, class PointIntp,DataHolder side=HOST>
	//void Test_FMM_Circle(real dx, Vector<real, d> center = Vector<real, d>::Zero(), real radius = (real)1) {
	//	Typedef_VectorD(d);
	//	Sphere<d> geom(center, radius);
	//	//First test the init and assign
	//	VectorDi counts = (Vector<real, d>::Ones() * radius * 2.2 / dx).template cast<int>();
	//	Grid<d> grid(counts, dx, center - Vector<real, d>::Ones() * radius * 1.1, MAC);
	//	Timer timer;
	//	timer.Reset();
	//	LevelSet<d, PointIntp, HOST> levelset_host;
	//	levelset_host.Init(grid, geom);
	//	timer.Record("init with geom");
	//	//Then I will set the none surface point to max
	//	LevelSet<d, PointIntp, HOST> levelset_tem(grid);
	//	levelset_tem.Init(levelset_host);
	//	Array<bool, HOST> is_surface(grid.DoF());
	//	for (int i = 0; i < levelset_tem.grid.DoF(); i++) {
	//		is_surface[i] = false;
	//		VectorDi cell = grid.Coord(i);
	//		for (int j = 0; j < 2 * d; j++) {
	//			VectorDi nb = grid.Neighbor_Node(cell, j);
	//			if (grid.Valid(nb) && levelset_tem.Interface(levelset_tem.Phi(grid.Position(cell)), levelset_tem.Phi(grid.Position(nb))) ) {
	//				is_surface[i] = true;
	//				break;
	//			}
	//		}
	//	}
	//	for (int i = 0; i < levelset_tem.grid.DoF(); i++) {
	//		if (!is_surface[i]) {
	//			(*(levelset_tem.phi.data))[i] = levelset_tem.Sign((*(levelset_tem.phi.data))[i])*std::numeric_limits<real>::max();
	//		}
	//	}
	//	LevelSet<d, PointIntp, side> levelset;
	//	levelset.Init(levelset_tem);
	//	timer.Record("break the distance function");
	//	//Here the bandwith=-1, with the meaning that there is no limit for distance
	//	Info("before fast marching phi: \n{}", levelset.phi);
	//	levelset.Fast_Marching(-1);
	//	timer.Record("fast marching");
	//	Info("after fast marching phi: \n{}", levelset.phi);
	//	levelset_tem.Init(levelset);
	//	real max_error = 0;
	//	for (int i = 0; i < levelset_tem.grid.DoF(); i++) {
	//		VectorDi cell = grid.Coord(i);
	//		VectorD pos = grid.Position(cell);

	//		bool continue_flag = false;
	//		for (int j = 0; j < d; j++) {
	//			if (cell[j] == 0 || cell[j] == grid.counts[j] - 1) continue_flag = true;
	//			if ((pos - center).norm() < 2 * dx) continue_flag = true;
	//		}
	//		if (continue_flag) continue;

	//		real refPhi = levelset_host.Phi(pos);
	//		real actPhi = levelset_tem.Phi(pos);
	//		//printf("(%f,%f)\t", refPhi, actPhi);
	//		Assert(fabs(refPhi - actPhi)< 2*dx, "Test_levelset_FMM failed: index {} failed,pos is {},{},ref phi is {},actual phi is {}", i, pos[0], pos[1], refPhi, actPhi);
	//		if (fabs(refPhi - actPhi)  > max_error) max_error = fabs(refPhi - actPhi);
	//	}
	//	timer.Record("compare phi one by one");
	//	//timer.Output_Profile();
	//	Pass("Test_FarstMarching_levelset_host (Circle) Passed, with max error {}!",max_error);

	//}
	//template<int d, class PointIntp>
	//void Test_FMM_Circle_GPU(real dx, Vector<real, d> center = Vector<real, d>::Zero(), real radius = (real)1) {

	//}
}