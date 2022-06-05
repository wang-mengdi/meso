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
		Field<bool, d> fixed;
		FaceField<T, d> vol;
		FaceField<bool, d> face_fixed;
		FaceField<T, d> initial_vel;

		void Apply(json& j, FluidFreeSurface<T, d>& fluid) {
			int test = Json::Value(j, "test", 0);
			switch (test) {
			case 0:Case_0(j, fluid); break;
			default:Assert(false, "test {} not exist", test); break;
			}
		}

		void Case_0(json& j, FluidFreeSurface<T, d>& fluid) {
			int scale = Json::Value<int>(j, "scale", 32);
			T side_len = 1.0;
			T dx = side_len / scale;
			VectorDi grid_size = scale * MathFunc::Vi<d>(1, 2, 1);
			Grid<d> grid(grid_size, dx, VectorD::Zero(), MAC);

			Eigen::Matrix<int, 3, 2> bc_width;
			Eigen::Matrix<real, 3, 2> bc_val;
			bc_width << 1, 1, 1, 1, 1, 1;
			bc_val << 0, 0, 0, 0, 0, 0;

			GridEulerFunc::Set_Boundary(grid, bc_width, bc_val, fixed, vol, face_fixed, initial_vel);

			Plane<d> plane(Vector<T, d>::Unit(1), grid.Center());

			fluid.Init(j, plane, fixed, vol, face_fixed, initial_vel);
		}
	};
}