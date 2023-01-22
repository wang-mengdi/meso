//////////////////////////////////////////////////////////////////////////
// Poisson Mapping Tests
// Copyright (c) (2022-), Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "DampedJacobiSmoother.h"
#include "PoissonFunc.h"
#include "Random.h"
#include "IOFunc.h"

namespace Meso {
	template<class T, int d>
	MaskedPoissonMapping<T, d> Random_Poisson_Mapping(const Grid<d> grid, real max_vol = 1) {
		Typedef_VectorD(d);
		FaceField<T, d> vol(grid);
		Field<unsigned char, d> cell_type(grid);
		MaskedPoissonMapping<T, d> mapping;
		vol.Iterate_Faces([&](const int axis, const VectorDi face) {vol(axis, face) = Random::Uniform(0, max_vol); });
		cell_type.Iterate_Nodes([&](const VectorDi cell) {	cell_type(cell) = (unsigned char)Random::RandInt(0, 2);	});
		mapping.Init(cell_type, vol);
		return mapping;
	}

	//vol==1, it's a 0/1 poisson mapping
	template<class T, int d>
	MaskedPoissonMapping<T, d> Random_Poisson01_Mapping(const Grid<d> grid) {
		Typedef_VectorD(d);
		FaceField<T, d> vol(grid);
		Field<unsigned char, d> cell_type(grid);
		MaskedPoissonMapping<T, d> mapping;
		vol.Iterate_Faces([&](const int axis, const VectorDi face) {vol(axis, face) = 1; });
		cell_type.Iterate_Nodes([&](const VectorDi cell) {	cell_type(cell) = (unsigned char)Random::RandInt(0, 1);	});
		mapping.Init(cell_type, vol);
		return mapping;
	}

	template<class T, int d>
	void Test_Poisson_Diagonal(Vector<int, d> counts) {
		Typedef_VectorD(d);
		typedef Eigen::Matrix<T, Eigen::Dynamic, 1> EigenVec;
		Grid<d> grid(counts);
		FaceField<T, d> vol(grid);
		Field<unsigned char, d> cell_type(grid);
		MaskedPoissonMapping<T, d> mapping;
		vol.Iterate_Faces(
			[&](const int axis, const VectorDi face) {
				vol(axis, face) = Random::Uniform(0, 1);
			}
		);
		cell_type.Iterate_Nodes(
			[&](const VectorDi cell) {
				cell_type(cell) = (unsigned char)Random::RandInt(0, 2);
			}
		);
		mapping.Init(cell_type, vol);
		ArrayDv<T> diag_dev(mapping.YDoF());
		PoissonLike_Diagonal(diag_dev, mapping);
		Array<T> diag_host = diag_dev;

		EigenVec vec_diag_grdt(mapping.YDoF());
		grid.Exec_Nodes(
			[&](const VectorDi cell) {
				int idx = grid.Index(cell);
		if (cell_type(cell) == 1 || cell_type(cell) == 2) {
			vec_diag_grdt[idx] = 1;
		}
		else {
			vec_diag_grdt[idx] = 0;
			for (int axis = 0; axis < d; axis++) {
				VectorDi face0 = cell, face1 = cell + VectorDi::Unit(axis);
				vec_diag_grdt[idx] += vol(axis, face0);
				vec_diag_grdt[idx] += vol(axis, face1);
			}
		}
			}
		);
		EigenVec vec_diag_host(mapping.YDoF());
		for (int i = 0; i < mapping.YDoF(); i++) vec_diag_host[i] = diag_host[i];
		if (vec_diag_grdt.isApprox(vec_diag_host)) Pass("Test_Poisson_Diagonal passed for counts={}", counts);
		else Error("Test_Poisson_Diagonal failed for counts={}", counts);
		//Info("vec_diag_grdt: {}", vec_diag_grdt);
		//Info("vec_diag_host: {}", vec_diag_host);
		//Info("error: {}", vec_diag_host - vec_diag_grdt);
	}


	template<class T, int d>
	void Test_Damped_Jacobian(int n) {
		Typedef_VectorD(d);
		Grid<d> grid(VectorDi::Ones() * n);
		n = grid.DoF();
		MaskedPoissonMapping<T, d> mapping = Random_Poisson_Mapping<T>(grid);
		ArrayDv<T> rhs = Random::Random_Array<T>(n, (T)0.0, (T)1.0);
		for (int i = 0; i < 100; i++) {
			DampedJacobiSmoother<T> smoother(mapping, i);
			ArrayDv<T> x(n), res(n);
			smoother.Apply(x, rhs);
			mapping.Apply(res, x);
			ArrayFunc::Minus(res, rhs);
			real l2_error = ArrayFunc::Norm(res);
			Info("iter {} l2_error {}", i, l2_error);
		}
	}

	template<class T, int d>
	void Test_Search_Boundary(const Vector<int, d> _counts, bool _output)
	{
		Typedef_VectorD(d);
		real domain_size = 1.0;
		Grid<d> grid(_counts, domain_size / _counts[0], -VectorD::Ones() * domain_size * 0.5);
		Field<unsigned char, d> cell_type(grid);
		grid.Iterate_Nodes([&](const VectorDi cell) {
			for (int axis = 0; axis < d; axis++)
			{
				if (cell[axis] == 0)
				{
					cell_type(cell) = 1;
					return;
				}
				else if (cell[axis] == grid.Counts()[axis] - 1)
				{
					cell_type(cell) = 2;
					return;
				}
			}
			cell_type(cell) = 0;
			});
		FaceField<real, d> vol(grid);
		MaskedPoissonMapping<T, d> poisson;
		poisson.Init(cell_type, vol);
		poisson.Search_Boundary();
		if (_output)
		{
			Field<unsigned char, d> cell_type_host = poisson.cell_type;
			std::string cell_type_name = "cell_type.vts";
			VTKFunc::Write_VTS(cell_type_host, cell_type_name);
			Array<int> boundary_tiles_host = poisson.boundary_tiles;
			for (auto a : boundary_tiles_host)
				Info("{}", a);
		}
	}

	template<class T, int d>
	void Test_Boundary_Apply_With_Apply(const Vector<int, d> _counts)
	{
		Typedef_VectorD(d);
		real domain_size = 1.0;
		Grid<d> grid(_counts, domain_size / _counts[0], -VectorD::Ones() * domain_size * 0.5);
		MaskedPoissonMapping<T, d> poisson = Random_Poisson_Mapping<T, d>(grid);
		poisson.Search_Boundary();
		Info("boundary tiles: {}", poisson.boundary_tiles.size());
		int n = grid.Memory_Size();
		Array<T> p_host(n);
		for (int i = 0; i < n; i++)
			p_host[i] = Random::Uniform(0, 1);
		ArrayDv<T> p_dev = p_host;
		ArrayDv<T> Ap_dev(n);
		ArrayDv<T> Ap_boundary_dev(n);
		poisson.Apply(Ap_dev, p_dev);
		poisson.Boundary_Apply(Ap_boundary_dev, p_dev);
		Array<T> Ap_host = Ap_dev;
		Array<T> Ap_boundary_host = Ap_boundary_dev;
		Field<unsigned char, d> cell_type_host = poisson.cell_type;
		grid.Iterate_Nodes([&](const VectorDi cell) {
			int cell_id = grid.Index(cell);
		if (cell_type_host(cell) == 3)
		{
			if (fabs(Ap_host[cell_id] - Ap_boundary_host[cell_id]) > 1e-10)
				Assert(false, "Poisson Boundary_Apply is not consistent with Apply");
		}
			});
			
	}
}