//////////////////////////////////////////////////////////////////////////
// Poisson Mapping Tests
// Copyright (c) (2022-), Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "PoissonSmoother.h"
#include "Random.h"

namespace Meso {

	template<class T, int d>
	void Test_Poisson_Diagonal(Vector<int, d> counts) {
		Typedef_VectorD(d);
		typedef Eigen::Matrix<T, Eigen::Dynamic, 1> EigenVec;
		Grid<d> grid(counts);
		FaceField<T, d> vol(grid);
		Field<bool, d> fixed(grid);
		PoissonMapping<T, d> mapping;
		vol.Iterate_Faces(
			[&](const int axis, const VectorDi face) {
				vol(axis, face) = Random::Uniform(0, 1);
			}
		);
		fixed.Iterate_Cells(
			[&](const VectorDi cell) {
				fixed(cell) = !(bool)Random::RandInt(0, 9);
			}
		);
		mapping.Init(grid, vol, fixed);
		ArrayDv<T> diag_dev(grid.DoF());
		Poisson_Diagonal(diag_dev, mapping);
		Array<T> diag_host = diag_dev;

		EigenVec vec_diag_grdt(grid.DoF());
		grid.Exec_Cells(
			[&](const VectorDi cell) {
				int idx = grid.Index(cell);
				if (fixed(cell)) {
					vec_diag_grdt[idx] = 1;
				}
				else {
					vec_diag_grdt[idx] = 0;
					for (int axis = 0; axis < d; axis++) {
						vec_diag_grdt[idx] += vol(axis, cell);
						vec_diag_grdt[idx] += vol(axis, cell + VectorDi::Unit(axis));
					}
				}
			}
		);
		EigenVec vec_diag_host(grid.DoF());
		for (int i = 0; i < grid.DoF(); i++) vec_diag_host[i] = diag_host[i];
		if (vec_diag_grdt.isApprox(vec_diag_host)) Pass("Test_Poisson_Diagonal passed for counts={}", counts);
		else Error("Test_Poisson_Diagonal failed for counts={}", counts);
		//Info("vec_diag_grdt: {}", vec_diag_grdt);
		//Info("vec_diag_host: {}", vec_diag_host);
		//Info("error: {}", vec_diag_host - vec_diag_grdt);
	}

}