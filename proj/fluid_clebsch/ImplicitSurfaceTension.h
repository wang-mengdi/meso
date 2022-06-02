//////////////////////////////////////////////////////////////////////////
// Implicitly solve surface tension
// Copyright (c) (2018-), Mengdi Wang
// This file is part of SimpleX,whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "MacGrid.h"
#include "FaceField.h"
#include "BoundaryCondition.h"
#include "LevelSet.h"
#include "Json.h"
//#include "SparseSolverCPX.h"
#include "Timer.h"
#include <fmt/os.h>

template<class T, int d>
class ImplicitSurfaceTension {
	Typedef_VectorDii(d);
public:
	MacGrid<d> mac_grid;
	const BoundaryConditionMacGrid<d>* bc = nullptr;
	const LevelSet<d>* levelset = nullptr;//everything in levelset is real type
	//FaceField<T, d> u_old;
	//Array<DiagonalPoissonSolver<Scalar, d>> solver_cpx;

	//FaceField<int, d> macgrid_to_matrix;
	//Array<std::pair<int, int>> matrix_to_macgrid;

	//Upon construction, mac_grid must be fully initialized
	ImplicitSurfaceTension(const Meso::json &j, const MacGrid<d> &_mac_grid, const BoundaryConditionMacGrid<d>* _bc, const LevelSet<d>* _levelset) {
		mac_grid = _mac_grid;
		bc = _bc;
		levelset = _levelset;
		//u_old.Resize(mac_grid.grid.cell_counts, 0);
		//solver_cpx.resize(d);
		//for (int axis = 0;axis < d;axis++) {
		//	VectorDi face_grid_counts = mac_grid.grid.cell_counts;
		//	face_grid_counts[axis] ++;
		//	MacGrid<d> face_mac_grid(face_grid_counts);
		//	solver_cpx[axis].Init(face_mac_grid, 100);
		//}
	}
	T Levelset_Dirac(const T phi, const T dirac_band_width) const
	{
		if (phi < -dirac_band_width) return 0;
		else if (phi > dirac_band_width) return 0;
		else return 0.5 * (1.0 + cos(pi * phi / dirac_band_width)) / dirac_band_width;
	}

	bool Is_Levelset_Interface_Face(const int axis, const VectorDi& face, const T narrow_band_width) const {
		if (!mac_grid.Valid_Face(axis, face) || bc->Is_Psi_N(axis, face)) return false;
		//check neighbor cell, filtering out most faces
		//NOTE: seems these lines will not filter much things
		//VectorDi nb_cell0 = MacGrid<d>::Face_Incident_Cell(axis, face, 0);
		//if (fabs(levelset->phi(nb_cell0)) > narrow_band_width + mac_grid.grid.dx * 0.5) return false;
		//VectorDi nb_cell1 = MacGrid<d>::Face_Incident_Cell(axis, face, 1);
		//if (fabs(levelset->phi(nb_cell1)) > narrow_band_width + mac_grid.grid.dx * 0.5) return false;
		//check the face itself
		VectorD pos = mac_grid.Face_Center(axis, face);
		auto phi = levelset->Phi(pos);
		return abs(phi) < narrow_band_width;
	}

	bool Is_Levelset_Interface_Face_Index(const std::pair<int, int> face_idx, const T narrow_band_width) const
	{
		int axis = face_idx.first;
		VectorDi face = mac_grid.face_grids[axis].Node_Coord(face_idx.second);
		VectorD pos = mac_grid.Face_Center(axis, face);
		if (bc->Is_Psi_N(axis, face)) return false;
		auto phi = levelset->Phi(pos);
		return abs(phi) < narrow_band_width;
	}

	////[NOTE] when I wrote the Solve() function, I thought face_grid[axis] holds all faces or direction axis as "cells"
	// However it seems that they're actually "nodes"
	// That may not cause too much trouble within a properly bounded system, but remember to check this
	// ----Mengdi
	
	////Sparse version, parallel build system
	void Solve(const T dt, FaceField<T, d>& velocity, const T sigma, const T narrow_band_width, const T dirac_band_width);

//	void Solve1(const T dt, FaceField<T, d>& velocity, const T sigma, const T narrow_band_width, const T dirac_band_width)
//	{
//		for (int axis = 0;axis < d;axis++) {
//			//for fixed term
//			auto is_unknown_func = std::bind(
//				&ImplicitSurfaceTension<T, d>::Is_Levelset_Interface_Face,
//				this,
//				axis,
//				std::placeholders::_1,
//				narrow_band_width
//			);
//			//for extra diagonal term
//			auto extra_diag_func = [&](const VectorDi& cell) {
//				return (Scalar)1.0;
//			};
//			//for laplacian term
//			//"face_face" means if you consider face nodes of axis as cells, the face between it
//			auto vol_func = [&](const int face_axis, const VectorDi& face_face)->T {
//				VectorD poses[2];
//				for (int i = 0;i < 2;i++) {
//					//the actual face node
//					VectorDi face = MacGrid<d>::Face_Incident_Cell(face_axis, face_face, i);
//					poses[i] = mac_grid.Face_Center(axis, face);
//				}
//				return sigma * Levelset_Dirac(levelset->Phi((poses[0] + poses[1]) * .5), dirac_band_width) * MathFunc::Power2(dt / mac_grid.grid.dx);
//			};
//
//			//fill u_old
//			int face_num = mac_grid.Number_Of_Faces(axis);
//#pragma omp parallel for
//			for (int it = 0; it < face_num; it++) {
//				VectorDi face = mac_grid.Face_Coord(axis, it);
//				if (!is_unknown_func(face)) {
//					u_old(axis, face) = 0;
//				}
//				else {
//					VectorD pos = mac_grid.Face_Center(axis, face);
//					auto kappa = levelset->Curvature(pos);
//					VectorD normal = levelset->Normal(pos);
//					u_old(axis, face) = velocity(axis, face) - sigma * Levelset_Dirac(levelset->Phi(pos), dirac_band_width) * dt * kappa * normal[axis];
//					//todo: these are not necessary, let 
//					for (int i = 0; i < Grid<d>::Number_Of_Nb_C(); i++) {
//						VectorDi nb_face = Grid<d>::Nb_C(face, i);
//						VectorD nb_pos = mac_grid.Face_Center(axis, nb_face);
//						if (mac_grid.Valid_Face(axis, nb_face) && !is_unknown_func(nb_face)) {
//							int face_axis;VectorDi face_face;
//							MacGrid<d>::Cell_Incident_Face(face, i, face_axis, face_face);
//							u_old(axis, face) += velocity(axis, nb_face) * vol_func(face_axis, face_face);
//						}
//					}
//				}
//			}
//
//			solver_cpx[axis].Update_System(is_unknown_func, extra_diag_func, vol_func);
//			solver_cpx[axis].Solve(u_old.face_fields[axis]);
//		}
//		for (int axis = 0;axis < d;axis++) {
//			int face_num = mac_grid.Number_Of_Faces(axis);
//#pragma omp parallel for
//			for (int it = 0; it < face_num; it++) {
//				VectorDi face = mac_grid.Face_Coord(axis, it);
//				PoissonDescriptor<d>& descr = solver_cpx[axis].diagonal_mapping.descr;
//				int cell_ind = -1;
//				if constexpr (d == 2) cell_ind = descr.grid.cell_ind(face[0], face[1]);
//				else if constexpr (d == 3) cell_ind = descr.grid.cell_ind(face[0], face[1], face[2]);
//				if (!descr.h_fixed[cell_ind]) {
//					velocity(axis, face) = solver_cpx[axis].x_field_host(face);
//				}
//			}
//		}
//	}
//	
//
//	////Pure sparse version, sequential build system
//	
//	void Solve2(const T dt, FaceField<T, d>& velocity, const T sigma, const T narrow_band_width, const T dirac_band_width)
//	{
//		Timer timer; timer.Reset();
//
//		// initialize fluid
//		std::function<bool(const std::pair<int, int>)> Is_Interface_Face_Index =
//			std::bind(
//				&ImplicitSurfaceTension<T, d>::Is_Levelset_Interface_Face_Index,
//				this,
//				std::placeholders::_1, narrow_band_width
//			);
//
//		matrix_to_macgrid.clear();
//		Build_MacGrid_Face_Matrix_Bijective_Mapping(
//			mac_grid, 
//			Is_Interface_Face_Index,
//			macgrid_to_matrix, 
//			matrix_to_macgrid
//		);
//		int n = matrix_to_macgrid.size();
//		if (n == 0) {
//			timer.Elapse_And_Output_And_Reset("Empty surface tension system, exit: ");
//			return;
//		}
//		SparseMatrixT B;
//		VectorX u_new;
//		VectorX u_old;
//
//		// setup A, x, and b
//		B.resize(n, n);
//		u_new.resize(n); u_new.fill(0);
//		u_old.resize(n); u_old.fill(0);
//		Array<TripletT> elements;
//
//		for (int r = 0; r < n; r++) {
//			const int axis = matrix_to_macgrid[r].first;
//			const VectorDi& face = mac_grid.face_grids[axis].Node_Coord(matrix_to_macgrid[r].second);
//			VectorD pos = mac_grid.Face_Center(axis, face);
//			auto kappa = levelset->Curvature(pos);
//			VectorD normal = levelset->Normal(pos);
//			T dia_coef = 1;
//			//neglect Jacobi and Hessian
//			u_old[r] = velocity(axis, face) - sigma * Levelset_Dirac(levelset->Phi(pos), dirac_band_width) * dt * kappa * normal[axis];
//			for (int i = 0; i < Grid<d>::Number_Of_Nb_C(); i++) {
//				VectorDi nb_face = Grid<d>::Nb_C(face, i);
//				VectorD nb_pos = mac_grid.Face_Center(axis, nb_face);
//				if (!mac_grid.Valid_Face(axis, nb_face) || bc->Is_Psi_N(axis, nb_face)) continue;
//				T a = sigma * Levelset_Dirac(levelset->Phi((pos + nb_pos) * .5), dirac_band_width) * dt * dt / (mac_grid.grid.dx * mac_grid.grid.dx);
//				dia_coef += a;
//				int c = macgrid_to_matrix(axis, nb_face);
//				if (Is_Interface_Face_Index(std::make_pair(axis, mac_grid.face_grids[axis].Node_Index(nb_face)))) {
//					elements.push_back(TripletT(r, c, -a));
//				}
//				else u_old[r] += a * velocity(axis, nb_face);
//			}
//			elements.push_back(TripletT(r, r, dia_coef));
//		}
//		B.setFromTriplets(elements.begin(), elements.end()); B.makeCompressed();
//		timer.Elapse_And_Output_And_Reset("Construct matrix for implicit surface tension: ");
//
//		SparseSolverCPX<Scalar> cpx_solver;
//		cpx_solver.Init(B, 3000, 1e-5);
//		cpx_solver.Solve(u_new.data(), u_old.data());
//		timer.Elapse_And_Output_And_Reset("Solve for implicit surface tension: ");
//		//#ifdef USE_CUDA
//		//MultiGridCuda::Preconditioned_Conjugate_Gradient<SparseMatrix<T>, T>(&B, &u_old[0], &u_new[0], 3000, (T)1e-5,/*diagonal precond*/true);	////GPU D-PCG
//		//#endif
//
//
//		#pragma omp parallel for
//		for (int r = 0; r < n; r++) {
//			const int axis = matrix_to_macgrid[r].first;
//			const VectorDi& face = mac_grid.face_grids[axis].Node_Coord(matrix_to_macgrid[r].second);
//			velocity(axis, face) = u_new[r];
//		}
//		timer.Elapse_And_Output_And_Reset("Correction for implicit surface tension: ");
//	}


};

