//////////////////////////////////////////////////////////////////////////
// Clebsch Fluid: Implementation of Schrödinger's Smoke [Chern et al. 2015]
// Copyright (c) (2022-), Bo Zhu, Mengdi Wang, Zhecheng Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once

#include <thrust/complex.h>
#include <array>

#include "Multigrid.h"
#include "ConjugateGradient.h"
#include "Advection.h"
#include "Simulator.h"
#include "IOFunc.h"
#include "device_launch_parameters.h"

namespace Meso {
	using namespace std::complex_literals;
	using Vector2C = Vector<thrust::complex<real>, 2>;
	__global__ void W2V_Mapping_Kernel2(const Grid<2> grid, real* face_x, real* face_y, const Vector2C* cell, const real h_bar_over_dx)
	{
		const int i = threadIdx.x + blockIdx.x * blockDim.x;
		const int j = threadIdx.y + blockIdx.y * blockDim.y;

		if (i > grid.counts[0] - 1 || j > grid.counts[1] - 1) return;

		const Vector2C cell_data = cell[grid.Index(i, j)];
		int face_ind;

		// x-axis faces
		face_ind = grid.Face_Index(0, Vector2i(i + 1, j));
		face_x[face_ind] = h_bar_over_dx * thrust::arg(cell_data.dot(cell[grid.Index(i + 1, j)]));

		// y-axis faces
		face_ind = grid.Face_Index(1, Vector2i(i, j + 1));
		face_y[face_ind] = h_bar_over_dx * thrust::arg(cell_data.dot(cell[grid.Index(i, j + 1)]));
	}

	template<int d>
	class FluidClebsch : public Simulator {
		Typedef_VectorD(d);
	public:
		real h_bar;
		FieldDv<Vector2C, d> wave_function;
		FaceFieldDv<real, d> velocity;

		BoundaryConditionDirect<FieldDv<Vector2C, d>> psi_D;
		FaceFieldDv<real, d> temp_velocity;
		FieldDv<real, d> pressure;
		FieldDv<real, d> vel_div;
		MaskedPoissonMapping<real, d> poisson;
		VCycleMultigridIntp<real, d> MG_precond;
		ConjugateGradient<real> MGPCG;
		void Init(real h_bar_, Field<bool, d>& fixed, Field<Vector2C, d>& initial_wave_function, FaceField<real, d>& vol) {
			h_bar = h_bar_;
			wave_function.Deep_Copy(initial_wave_function);
			Wave_Function_To_Velocity(velocity, wave_function);
			psi_D.Init(fixed, initial_wave_function);
			temp_velocity.Init(velocity.grid);
			pressure.Init(velocity.grid);
			vel_div.Init(velocity.grid);

			poisson.Init(fixed, vol);
			MG_precond.Init_Poisson(poisson, 2, 2);
			MGPCG.Init(&poisson, &MG_precond, false, -1, 1e-6);
		}
		virtual real CFL_Time(const real cfl) {
			real dx = velocity.grid.dx;
			real max_vel = velocity.Max_Abs();
			return dx * cfl / max_vel;
		}
		virtual void Output(DriverMetaData& metadata) {
			std::string vts_name = fmt::format("vts{:04d}.vts", metadata.current_frame);
			bf::path vtk_path = metadata.base_path / bf::path(vts_name);
			VTKFunc::Write_VTS(velocity, vtk_path.string());
		}
		virtual void Wave_Function_To_Velocity(FaceFieldDv<real, d>& F, const FieldDv<Vector2C, d>& C) {
			const real h_bar_over_dx = h_bar / C.grid.dx;
			F.Init(C.grid, MathFunc::Zero<real>());
			dim3 blocknum, blocksize;
			C.grid.Get_Kernel_Dims(blocknum, blocksize);
			const Vector2C* cell = C.Data_Ptr();
			if constexpr (d == 2) {
				real* face_x = F.Data_Ptr(0);
				real* face_y = F.Data_Ptr(1);

				W2V_Mapping_Kernel2 << <blocknum, blocksize >> > (C.grid, face_x, face_y, cell, h_bar_over_dx);
			}
			/*else if constexpr (d == 3) {
				real* face_x = F.Data_Ptr(0);
				real* face_y = F.Data_Ptr(1);
				real* face_z = F.Data_Ptr(2);

				W2V_Mapping_Kernel3 << <blocknum, blocksize >> > (C.grid, face_x, face_y, face_z, cell, h_bar_over_dx);
			}*/
		}
		virtual void Advance(DriverMetaData& metadata) {
			real dt = metadata.dt;

			//TODO: solve
			

			//wave function to velocity
			Wave_Function_To_Velocity(velocity, wave_function);

			//projection
			//vel_div=div(velocity)
			Exterior_Derivative(vel_div, velocity);

			auto [iter, res] = MGPCG.Solve(pressure.Data(), vel_div.Data());
			Info("Solve poisson with {} iters and residual {}", iter, res);

			//velocity+=grad(p)
			Exterior_Derivative(temp_velocity, pressure);
			temp_velocity *= poisson.vol;

			velocity += temp_velocity;
			psi_D.Apply(wave_function);

			Exterior_Derivative(vel_div, velocity);
			Info("After projection max div {}", vel_div.Max_Abs());
		}
	};
}