//////////////////////////////////////////////////////////////////////////
// Level set
// Copyright (c) (2018-), Bo Zhu, Xingyu Ni, Fan Feng
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#include <numeric>
#include <set>
#include <queue>
#include <utility>
#include <iostream>
#include "LevelSet.h"
#include "Constants.h"
#include "Timer.h"

namespace Meso {

	real Solve_Upwind_Eikonal2(const real p1, const real p2, const real dx)
	{
		if (abs(p1) >= abs(p2) + dx)  return p2 + dx;
		else if (abs(p2) >= abs(p1) + dx)  return p1 + dx;
		else {
			real delta = (real)2 * dx * dx - pow(p1 - p2, 2);
			Assert(delta >= 0, "Error: solve quadratic delta={}", delta);
			return (real).5 * (p1 + p2 + sqrt(delta));
		}
	}

	real Solve_Upwind_Eikonal3(const real p1, const real p2, const real p3, const real dx)
	{
		real delta = pow(p1 + p2 + p3, 2) - (real)3 * (p1 * p1 + p2 * p2 + p3 * p3 - dx * dx);
		if (delta < (real)0) {
			int i = 0; real p_max = abs(p1); if (abs(p2) > p_max) { i = 1; p_max = abs(p2); }if (abs(p3) > p_max) { i = 2; p_max = abs(p3); }
			real q1, q2; if (i == 0) { q1 = p2; q2 = p3; }
			else if (i == 1) { q1 = p1; q2 = p3; }
			else { q1 = p1; q2 = p2; }
			return Solve_Upwind_Eikonal2(q1, q2, dx);
		}
		return (p1 + p2 + p3 + sqrt(delta)) / 3.0;
	}

	template<int d> LevelSet<d>::LevelSet(const Grid<d> _grid)
	{
		Init(_grid);
	}

	template<int d> real LevelSet<d>::Cell_Fraction(const VectorDi& cell) const
	{
		real dx = phi.grid.dx;
		return (real).5 - MathFunc::Clamp(phi(cell), -(real).5 * dx, (real).5 * dx) / dx;
	}

	//////////////////////////////////////////////////////////////////////////
	////Fast marching method

	template<int d> void LevelSet<d>::Fast_Marching(const real band_width, bool reinit_interface)
	{
		Grid<d> grid = phi.grid;		
		Field<real, d> tent(grid, band_width < 0 ? std::numeric_limits<real>::max() : band_width);
		Field<bool, d> done(grid, 0);	//Fan: can change it to a field
		
		std::priority_queue<PRI, Array<PRI>, std::greater<PRI> > heaps[2];
		const int cell_num = grid.Counts().prod();
		//real far_from_intf_phi_val=grid.dx*(real)5;

		//// Step 1: find interface cells and calculate the rough distance to interface
		grid.Exec_Nodes(
			[&](const VectorDi cell) {
				int idx=grid.Index(cell);
				for (int i = 0; i < Grid<d>::Neighbor_Node_Number(); i++) {
					VectorDi nb = grid.Neighbor_Node(cell, i);
					if (!grid.Valid(nb))continue;
					if (Is_Interface(cell, nb)) {
						if(reinit_interface){
							VectorD correct_phi = VectorD::Ones() * std::numeric_limits<real>::max();
							VectorDi correct_axis = VectorDi::Zero();
							for (int j = 0; j < Grid<d>::Neighbor_Node_Number(); j++) {
								VectorDi nb = grid.Neighbor_Node(cell, j);
								if (!grid.Valid(nb)) continue;
								const int nb_idx = grid.Index(nb);
								if (Is_Interface(cell, nb)) {
									real c_phi = Theta(phi(cell), phi(nb)) * grid.dx; // always non-negative
									int axis = grid.Neighbor_Node_Axis(j);
									correct_axis[axis] = 1;
									correct_phi[axis] = std::min(correct_phi[axis], c_phi);
								}
							}
							if (correct_axis != VectorDi::Zero()) {
								real hmnc_mean = (real)0;
								for (int k = 0; k < d; k++) {
									if (correct_axis[k] == 0)continue;
									hmnc_mean += (real)1 / (correct_phi[k] * correct_phi[k]);
								}
								hmnc_mean = sqrt((real)1 / hmnc_mean);
								tent(cell) = hmnc_mean;
							}
							else {
								Error("[Levelset] bad preconditioning");
							}
						}
						else {
							tent(cell) = abs(phi(cell));
						}
	#pragma omp critical
						{heaps[MathFunc::Sign(phi(cell)) > 0 ? 0 : 1].push(PRI(tent(cell), idx)); }
						
						break;
					}
				}
			}
		);
		
		//// Step 2: relax the field
#pragma omp parallel for
		for (int h = 0; h < 2; h++) {
			Relax_Heap(heaps[h], tent, done, phi);
		}

		ArrayFunc::Binary_Transform(
			phi.Data(),
			tent.Data(),
			[=](const real phi_i, const real tent_i) {return MathFunc::Sign(phi_i) * tent_i; },
			phi.Data()
		);
	}

	template class LevelSet<2>;
	template class LevelSet<3>;

}