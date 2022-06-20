//////////////////////////////////////////////////////////////////////////
// Initializer of Fluid Euler with Free Surface
// Copyright (c) (2022-), Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "FluidFreeSurface.h"
#include "ImplicitGeometry.h"
#include "Json.h"
#include "GridEulerFunc.h"

namespace Meso {
	template<class T, int d>
	class FluidFreeSurfaceInitializer {
	public:
		Typedef_VectorD(d);

		//define the boundary conditions of the eulerian system
		Field<CellType, d> cell_type;
		FaceField<T, d> initial_vel;

		void Apply(json& j, FluidFreeSurface<T, d>& fluid) {
			std::string test = Json::Value(j, "test", std::string("hydrostatic"));
			if (test == "hydrostatic") Init_Hydrostatic(j, fluid);
			else if (test == "dropfall") Init_DropFall(j, fluid);
			else Assert(false, "test {} not exist", test);
		}

		void Init_Hydrostatic(json& j, FluidFreeSurface<T, d>& fluid) {
			int scale = Json::Value<int>(j, "scale", 32);
			T side_len = 1.0;
			T dx = side_len / scale;
			VectorDi grid_size = scale * MathFunc::Vi<d>(1, 2, 1);
			Grid<d> grid(grid_size, dx, VectorD::Zero(), CENTER);

			cell_type.Init(grid, FLUID);
			Eigen::Matrix<int, 3, 2> bc_width;
			bc_width << 1, 1, 1, 1, 1, 1;
			GridEulerFunc::Set_Boundary_Cells(cell_type, bc_width, SOLID);
			initial_vel.Init(grid, 0);

			Plane<d> plane(Vector<T, d>::Unit(1), grid.Center());

			fluid.Init(j, plane, cell_type, initial_vel);
		}

		void Init_DropFall(json& j, FluidFreeSurface<T, d>& fluid) {
			int scale = Json::Value<int>(j, "scale", 32);
			T side_len = 1.0;
			T dx = side_len / scale;
			VectorDi grid_size = scale * MathFunc::Vi<d>(1, 2, 1);
			Grid<d> grid(grid_size, dx, VectorD::Zero(), CENTER);

			cell_type.Init(grid, FLUID);
			Eigen::Matrix<int, 3, 2> bc_width;
			bc_width << 1, 1, 1, 1, 1, 1;
			GridEulerFunc::Set_Boundary_Cells(cell_type, bc_width, SOLID);
			initial_vel.Init(grid, 0);

			Sphere<d> sphere(grid.Center() + VectorD::Unit(1) * 0.25, side_len * 0.2);

			fluid.Init(j, sphere, cell_type, initial_vel);
		}
	};
}