//////////////////////////////////////////////////////////////////////////
// Level set
// Copyright (c) (2018-), Bo Zhu, Xingyu Ni
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

	template<int d> void LevelSet<d>::Init(const Grid<d> _grid)
	{
		phi.Init(_grid, std::numeric_limits<real>::max());
	}


	//template<int d> real LevelSet<d>::Phi(const VectorD& pos) const
	//{
	//	return intp->Interpolate_Centers(phi, pos);
	//}

	//template<int d> real LevelSet<d>::Curvature(const VectorD& pos) const
	//{
	//	real one_over_dx = (real)1 / grid.dx; real one_over_two_dx = (real).5 * one_over_dx; real curvature = (real)0;
	//	for (int i = 0; i < d; i++) {
	//		VectorD normal_left = Normal(pos - VectorD::Unit(i) * grid.dx);
	//		VectorD normal_right = Normal(pos + VectorD::Unit(i) * grid.dx);
	//		curvature += (normal_right[i] - normal_left[i]) * one_over_two_dx;
	//	}
	//	return abs(curvature) < one_over_dx ? curvature : (curvature <= (real)0 ? (real)-1 : (real)1) * one_over_dx;

	//}

	//template<int d> Vector<real, d> LevelSet<d>::Closest_Point(const VectorD& pos, const real epsilon) const
	//{
	//	VectorD normal = Gradient(pos); normal.normalize();
	//	return pos - normal * (Phi(pos) + epsilon);
	//}

	//template<int d> Vector<real, d> LevelSet<d>::Closest_Point_With_Iterations(const VectorD& pos, const int max_iter/*=5*/) const
	//{
	//	VectorD intf_pos = pos;
	//	for (int i = 0; i < max_iter; i++) {
	//		intf_pos = Closest_Point(intf_pos);
	//		if (Phi(intf_pos) < (real)0)return intf_pos;
	//	}
	//	return intf_pos;
	//}

	template<int d> real LevelSet<d>::Cell_Fraction(const VectorDi& cell) const
	{
		real dx = phi.grid.dx;
		return (real).5 - MathFunc::Clamp(phi(cell), -(real).5 * dx, (real).5 * dx) / dx;
	}

	//////////////////////////////////////////////////////////////////////////
	////Fast marching method

	template<int d> void LevelSet<d>::Fast_Marching(const real band_width)
	{
		Grid<d> grid = phi.grid;
		//Timer timer;
		//timer.Reset();
		

		Field<real, d> tent(grid, band_width < 0 ? std::numeric_limits<real>::max() : band_width);
		Array<ushort> done(grid.DoF(), 0);
		using PRI = std::pair<real, int>;
		std::priority_queue<PRI, Array<PRI>, std::greater<PRI> > heaps[2];
		const int cell_num = grid.DoF();
		//real far_from_intf_phi_val=grid.dx*(real)5;

		//// Step 1: find interface cells
#pragma omp parallel for
		for (int i = 0; i < cell_num; i++) {
			const VectorDi cell = grid.Coord(i);
			//if(abs(phi(cell))>far_from_intf_phi_val)continue;		////ATTENTION: this might cause problem if the levelset is badly initialized

			for (int j = 0; j < Grid<d>::Neighbor_Node_Number(); j++) {
				VectorDi nb = grid.Neighbor_Node(cell, j);
				if (!grid.Valid(nb))continue;
				if (Is_Interface(cell, nb)) {
					done[i] = true; break;
				}
			}
		}
		//if (verbose)timer.Elapse_And_Output_And_Reset("FMM Precond: find interface");

		//// Step 2: calculate initial phi values for interface cells
#pragma omp parallel for
		for (int c = 0; c < cell_num; c++) {
			if (!done[c])continue;		////select interface cells
			const VectorDi cell = grid.Coord(c);

			VectorD correct_phi = VectorD::Ones() * std::numeric_limits<real>::max();
			VectorDi correct_axis = VectorDi::Zero();
			for (int i = 0; i < Grid<d>::Neighbor_Node_Number(); i++) {
				VectorDi nb = grid.Neighbor_Node(cell, i);
				if (!grid.Valid(nb)) continue;
				const int nb_idx = grid.Index(nb);
				if (done[nb_idx] && Is_Interface(cell, nb)) {
					real c_phi = Theta(phi(cell), phi(nb)) * grid.dx; // always non-negative
					int axis = grid.Neighbor_Node_Axis(i);
					correct_axis[axis] = 1;
					correct_phi[axis] = std::min(correct_phi[axis], c_phi);
				}
			}
			if (correct_axis != VectorDi::Zero()) {
				real hmnc_mean = (real)0;
				for (int i = 0; i < d; i++) {
					if (correct_axis[i] == 0)continue;
					hmnc_mean += (real)1 / (correct_phi[i] * correct_phi[i]);
				}
				hmnc_mean = sqrt((real)1 / hmnc_mean);
				tent(cell) = hmnc_mean;
			}
			else {
				Error("[Levelset] bad preconditioning");
			}
		}

		//// Step 3: perform relaxation on interface cells to fix their values

		//// initialize heap with front cells
#pragma omp parallel for
		for (int i = 0; i < cell_num; i++) {
			const VectorDi cell = grid.Coord(i);
			if (!done[i]) continue;
			auto [relax_success, val] = Relax_Node(cell, phi, tent, done);
#pragma omp critical
			{heaps[MathFunc::Sign(phi(cell)) > 0 ? 0 : 1].push(PRI(tent(cell), i)); }
		}

		//if (verbose)timer.Elapse_And_Output_And_Reset("FMM: Build heap");

		//// Step 4: relax the other part of field
#pragma omp parallel for
		for (int h = 0; h < 2; h++) {
			auto& heap = heaps[h];
			while (!heap.empty()) {
				const real top_val = heap.top().first;
				const int cell_idx = heap.top().second;
				const VectorDi cell = grid.Coord(cell_idx);
				heap.pop();
				if (tent(cell) != top_val) continue;
				done[cell_idx] = true;

				for (int i = 0; i < Grid<d>::Neighbor_Node_Number(); i++) {
					VectorDi nb = grid.Neighbor_Node(cell, i);
					if (!grid.Valid(nb))continue;
					const int nb_idx = grid.Index(nb);
					//relaxation
					if (!done[nb_idx]) {
						auto [relaxed, val] = Relax_Node(nb, phi, tent, done);
						if (relaxed) heap.push(PRI(val, nb_idx));
					}
				}
			}
		}

		//if (verbose)timer.Elapse_And_Output_And_Reset("FMM: Traverse heap");

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