//////////////////////////////////////////////////////////////////////////
// Level set
// Copyright (c) (2022-), Zhiqi Li, Bo Zhu, Xingyu Ni, Fan Feng
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "Grid.h"
#include "Field.h"
#include "Common.h"
#include "AuxFunc.h"
#include "Interpolation.h"
#include "ImplicitManifold.h"
#include<iostream>
//#include "cuda_runtime.h"
#include <typeinfo>
/// Some function is called as a intergral and some function is called only for a element
namespace Meso {
	real Solve_Upwind_Eikonal2(const real p1, const real p2, const real dx);
	real Solve_Upwind_Eikonal3(const real p1, const real p2, const real p3, const real dx);

	template<int d> 
	class LevelSet
	{
		Typedef_VectorD(d);
		using PRI = std::pair<real, int>;
	public:
		Field<real, d> phi;
		bool verbose = true;
	public:
		LevelSet() {}
		LevelSet(const Grid<d> _grid);
		void Init(const Grid<d> _grid) {
			phi.Init(_grid, std::numeric_limits<real>::max());
		}
		template<class Shape>
		void Init(const Grid<d> _grid, Shape geom) {
			Init(_grid);
			phi.Calc_Nodes(
				[&](const VectorDi cell) {
					return geom.Phi(_grid.Position(cell));
				}
			);
		}
		void Init(const Grid<d> _grid, const ImplicitManifold<d> &geom) {
			Init(_grid);
			phi.Calc_Nodes(
				[&](const VectorDi cell) {
					return geom.Phi(_grid.Position(cell));
				}
			);
		}

		void Init(const Field<real, d>& _phi) {
			phi=_phi;
		}

		real Phi(const VectorD pos) const {
			return IntpLinearClamp::Value(phi, pos);
		}

		real Cell_Fraction(const VectorDi& cell) const;		////approximate cell volume using phi value

		////Helper functions
		bool Is_Interface(const VectorDi cell0, const VectorDi cell1) const { return MathFunc::Sign(phi(cell0)) != MathFunc::Sign(phi(cell1)); }
		static real Theta(const real phi_1, const real phi_2) { return phi_1 / (phi_1 - phi_2); }
		
		//At the beginning of fast marching, the signs must be correct
		void Fast_Marching(const real band_width = (real)-1, bool reinit_interface = false);
	protected:
		//return [is_relaxed, cell_value]
		std::tuple<bool, real> Relax_Node(const VectorDi& cell, const Field<real, d>& phi, Field<real, d>& tent, const Field<bool,d>& done)
		{
			const Grid<d> grid = phi.grid;

			int sgn = MathFunc::Sign(phi(cell));

			// calculate correct phi from nb interface cells
			VectorD correct_phi = VectorD::Ones() * std::numeric_limits<real>::max();
			VectorDi correct_axis = VectorDi::Zero();
			for (int i = 0; i < Grid<d>::Neighbor_Node_Number(); i++) {
				VectorDi nb = grid.Neighbor_Node(cell, i);
				if (!grid.Valid(nb)) continue;
				if (done(nb)) {
					int axis = grid.Neighbor_Node_Axis(i); correct_axis[axis] = 1;
					correct_phi[axis] = std::min(correct_phi[axis], tent(nb));
				}
			}
			// update phi on the cell
			real new_phi;
			int n = correct_axis.sum();

			real old_tent = tent(cell);

			switch (n) {
			case 0: return std::make_tuple(false, old_tent);
			case 1: {
				real c_phi;
				for (int i = 0; i < d; i++)
					if (correct_axis[i] != 0) { c_phi = correct_phi[i]; break; }
				new_phi = grid.dx + c_phi;
			} break;
			case 2: {
				real p[2];
				int j = 0;
				for (int i = 0; i < d; i++)
					if (correct_axis[i] != 0) p[j++] = correct_phi[i];
				new_phi = Solve_Upwind_Eikonal2(p[0], p[1], grid.dx);
			} break;
			case 3: {
				new_phi = Solve_Upwind_Eikonal3(correct_phi[0], correct_phi[1], correct_phi[2], grid.dx);
			} break;
			default: {
				Error("[Levelset] bad solving Eikonal");
			} break;
			}

			
			if (new_phi < old_tent) {
				tent(cell) = new_phi;
				return std::make_tuple(true, new_phi);
			}
			else return std::make_tuple(false, old_tent);
		}

		//a normal Dijkstra
		void Relax_Heap(std::priority_queue<PRI, Array<PRI>, std::greater<PRI> > &heap, Field<real,d> &tent, Field<bool, d>& done, const Field<real,d> phi) {
			const Grid<d> grid = phi.grid;
			while (!heap.empty()) {
				const real top_val = heap.top().first;
				const int cell_idx = heap.top().second;
				const VectorDi cell = grid.Coord(cell_idx);
				heap.pop();
				if (tent(cell) != top_val) continue;
				done(cell) = true;

				for (int i = 0; i < Grid<d>::Neighbor_Node_Number(); i++) {
					VectorDi nb = grid.Neighbor_Node(cell, i);
					if (!grid.Valid(nb))continue;
					const int nb_idx = grid.Index(nb);
					if (done(nb)) continue;
					if (MathFunc::Sign(phi(cell)) != MathFunc::Sign(phi(nb))) continue;
					auto [relaxed, val] = Relax_Node(nb, phi, tent, done);
					if (relaxed) heap.push(PRI(val, nb_idx));
				}
			}
		}
	};
}