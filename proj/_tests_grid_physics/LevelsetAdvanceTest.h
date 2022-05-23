//////////////////////////////////////////////////////////////////////////
// Test leveset on GPU and this is advanced test by senario examples
// Copyright (c) (2022-), Zhiqi Li
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "LevelSet.h"
#include "Timer.h"

namespace Meso {
	/// Here, we test the fast marching method
	template<int d, class PointIntp,DataHolder side=HOST>
	void Test_FMM_Circle(real dx, Vector<real, d> center = Vector<real, d>::Zero(), real radius = (real)1) {
		Typedef_VectorD(d);
		Sphere<d> geom(center, radius);
		//First test the init and assign
		VectorDi counts = (Vector<real, d>::Ones() * radius * 2.2 / dx).template cast<int>();
		Grid<d> grid(counts, dx, center - Vector<real, d>::Ones() * radius * 1.1, MAC);
		Timer timer;
		timer.Reset();
		LevelSet<d, PointIntp, HOST> levelset_host(grid);
		levelset_host.Set_By_Geom(geom);
		timer.Record("init with geom");
		//Then I will set the none surface point to max
		LevelSet<d, PointIntp, HOST> levelset_tem(grid);
		levelset_tem.Initialize(levelset_host);
		Array<bool, HOST> is_surface(grid.DoF());
		for (int i = 0; i < levelset_tem.grid.DoF(); i++) {
			is_surface[i] = false;
			VectorDi cell = grid.Coord(i);
			for (int j = 0; j < 2 * d; j++) {
				VectorDi nb = grid.Nb_C(cell, j);
				if (grid.Valid(nb) && levelset_tem.Interface(levelset_tem.Phi(grid.Position(cell)), levelset_tem.Phi(grid.Position(nb))) ) {
					is_surface[i] = true;
					break;
				}
			}
		}
		for (int i = 0; i < levelset_tem.grid.DoF(); i++) {
			if (!is_surface[i]) {
				(*(levelset_tem.phi.data))[i] = levelset_tem.Sign((*(levelset_tem.phi.data))[i])*std::numeric_limits<real>::max();
			}
		}
		LevelSet<d, PointIntp, side> levelset;
		levelset.Initialize(levelset_tem);
		timer.Record("break the distance function");
		//Here the bandwith=-1, with the meaning that there is no limit for distance
		levelset.Fast_Marching(-1);
		timer.Record("fast marching");
		levelset_tem.Initialize(levelset);
		real max_error = 0;
		for (int i = 0; i < levelset_tem.grid.DoF(); i++) {
			VectorDi cell = grid.Coord(i);
			VectorD pos = grid.Position(cell);

			bool continue_flag = false;
			for (int j = 0; j < d; j++) {
				if (cell[j] == 0 || cell[j] == grid.counts[j] - 1) continue_flag = true;
				if ((pos - center).norm() < 2 * dx) continue_flag = true;
			}
			if (continue_flag) continue;

			real refPhi = levelset_host.Phi(pos);
			real actPhi = levelset_tem.Phi(pos);
			//printf("(%f,%f)\t", refPhi, actPhi);
			Assert(fabs(refPhi - actPhi)< 2*dx, "Test_levelset_FMM failed: index {} failed,pos is {},{},ref phi is {},actual phi is {}", i, pos[0], pos[1], refPhi, actPhi);
			if (fabs(refPhi - actPhi)  > max_error) max_error = fabs(refPhi - actPhi);
		}
		timer.Record("compare phi one by one");
		timer.Output_Profile();
		Pass("Test_FarstMarching_levelset_host (Circle) Passed, with max error {}!",max_error);

	}
	template<int d, class PointIntp>
	void Test_FMM_Circle_GPU(real dx, Vector<real, d> center = Vector<real, d>::Zero(), real radius = (real)1) {

	}
}