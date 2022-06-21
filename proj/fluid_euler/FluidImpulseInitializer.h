#pragma once

#include "FluidImpulse.h"
#include "ImplicitManifold.h"
#include "Json.h"

namespace Meso {

	template<int d>
	class  FluidImpulseInitializer {
	public:
		Typedef_VectorD(d);

		//define the boundary conditions of the eulerian system
		Field<bool, d> fixed;
		FaceField<real, d> vol;
		FaceField<bool, d> face_fixed;
		FaceField<real, d> initial_vel;

		void Apply(json& j, FluidImpulse<d>& fluid) {
			int test = Json::Value(j, "test", 0);
			switch (test) {
			case 0:Case_0(j, fluid); break;
			case 1:Case_1(j, fluid); break;
			case 2:Case_2(j, fluid); break;
			default:Assert(false, "test {} not exist", test); break;
			}
		}

		//an empty field with inflow boundary condition v=1
		void Case_0(json& j, FluidImpulse<d>& fluid) {
			int scale = Json::Value(j, "scale", 32);
			real dx = 1.0 / scale;
			VectorDi grid_size = scale * MathFunc::Vi<d>(2, 1, 1);
			Grid<d> grid(grid_size, dx, VectorD::Zero(), CENTER);

			Eigen::Matrix<int, 3, 2> bc_width;
			Eigen::Matrix<real, 3, 2> bc_val;
			bc_width << 0, -1, 0, 0, 0, 0;
			bc_val << 1, 1, 0, 0, 0, 0;

			GridEulerFunc::Set_Boundary(grid, bc_width, bc_val, fixed, vol, face_fixed, initial_vel);
			fluid.Init(fixed, vol, face_fixed, initial_vel);
			ArrayFunc::Fill(fluid.velocity.Data(0), 1.0);
		}

		void Case_1(json& j, FluidImpulse<d>& fluid) {
			int scale = Json::Value(j, "scale", 32);
			real dx = 1.0 / scale;
			VectorDi grid_size = scale * MathFunc::Vi<d>(1, 1, 1);
			Grid<d> grid(grid_size, dx, VectorD::Zero(), CENTER);

			Eigen::Matrix<int, 3, 2> bc_width;
			Eigen::Matrix<real, 3, 2> bc_val;
			bc_width << 1, 1, 1, 1, 1, 1;
			bc_val << 0, 0, 0, 0, 0, 0;

			GridEulerFunc::Set_Boundary(grid, bc_width, bc_val, fixed, vol, face_fixed, initial_vel);

			grid.Exec_Faces(
				[&](const int axis, const VectorDi face) {
					if (axis == 0 && face[1] >= grid.Counts()[1] - 2) {
						if (!face_fixed(axis, face)) {
							face_fixed(axis, face) = true;
							initial_vel(axis, face) = 1.0;
						}
					}
				}
			);

			fluid.Init(fixed, vol, face_fixed, initial_vel);
		}

		// with a sheprical obstacle, initial velocity of 1
		void Case_2(json& j, FluidImpulse<d>& fluid) {
			int scale = Json::Value(j, "scale", 32);
			real dx = 1.0 / scale;
			VectorDi grid_size = scale * MathFunc::Vi<d>(2, 1, 1);
			Grid<d> grid(grid_size, dx, VectorD::Zero(), CENTER);

			Eigen::Matrix<int, 3, 2> bc_width;
			Eigen::Matrix<real, 3, 2> bc_val;
			bc_width << 0, 0, 0, 0, 0, 0;
			bc_val << 1, 1, 0, 0, 0, 0;

			////sphere
			VectorD center = grid.Center();
			center[0] = center[0] / 2;
			real r = (real)0.1;
			Sphere<d> sphere(center, r);

			GridEulerFunc::Set_Boundary(grid, bc_width, bc_val, fixed, vol, face_fixed, initial_vel);

			grid.Exec_Nodes(
				[&](const VectorDi cell) {
					const VectorD pos = grid.Position(cell);
					if (sphere.Inside(pos)) {
						fixed(cell) = true;
					}
				}
			);

			grid.Exec_Faces(
				[&](const int axis, const VectorDi face) {
					const VectorD pos = grid.Face_Center(axis, face);
					if (sphere.Inside(pos)) {
						face_fixed(axis, face) = true;
						initial_vel(axis, face) = 0.0;
						vol(axis, face) = 0;
					}
				}
			);

			fluid.Init(fixed, vol, face_fixed, initial_vel);
		}
	};

}