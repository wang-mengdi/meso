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
#include "GridEulerFunc.h"
#include "device_launch_parameters.h"

namespace Meso {
	using namespace std::complex_literals;
	using Vector2C = Vector<thrust::complex<real>, 2>;
	__global__ void W2V_Mapping_Kernel2_Padding0(const Grid<2> grid, real* face_x, real* face_y, const Vector2C* cell, const real h_bar_over_dx)
	{
		Typedef_VectorD(2);
		VectorDi coord = GPUFunc::Thread_Coord<2>(blockIdx, threadIdx);

		if (!grid.Valid(coord)) return;

		const Vector2C cell_data = cell[grid.Index(coord)];
		const Vector2C cell_data_conj(thrust::conj(cell_data[0]), thrust::conj(cell_data[1]));

		int face_ind_x = Neighbor_Face_Index(grid, coord, 0, 1);
		int face_ind_y = Neighbor_Face_Index(grid, coord, 1, 1);

		if (face_ind_x != -1 && grid.Valid(coord + VectorDi::Unit(0))) {
			const Vector2C nb_cell_data = cell[grid.Index(coord + VectorDi::Unit(0))];
			face_x[face_ind_x] = h_bar_over_dx * thrust::arg(cell_data_conj.dot(nb_cell_data));
		}
		if (face_ind_y != -1 && grid.Valid(coord + VectorDi::Unit(1))) {
			const Vector2C nb_cell_data = cell[grid.Index(coord + VectorDi::Unit(1))];
			face_x[face_ind_y] = h_bar_over_dx * thrust::arg(cell_data_conj.dot(nb_cell_data));
		}
	}

	template<int d>
	__global__ void Wave_Function_Correction_Kernel(const Grid<d> grid, Vector2C* wave_function, const real* pressure, const real dx_over_h_har) {
		const int index = grid.Index(GPUFunc::Thread_Coord<d>(blockIdx, threadIdx));
		const real cell_p = pressure[index];
		Vector2C psi = wave_function[index];
		thrust::complex<real> c(thrust::exp(thrust::complex<real>(0, -1) * cell_p * dx_over_h_har));
		psi[0] *= c;
		psi[1] *= c;
		wave_function[index] = psi;
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
			Exterior_Derivative_W2V(velocity, wave_function);
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
			real max_vel = GridEulerFunc::Linf_Norm(velocity);
			return dx * cfl / max_vel;
		}
		virtual void Output(DriverMetaData& metadata) {
			std::string vts_name = fmt::format("vts{:04d}.vts", metadata.current_frame);
			bf::path vtk_path = metadata.base_path / bf::path(vts_name);
			VTKFunc::Write_VTS(velocity, vtk_path.string());
		}
		virtual void Schrödinger_Solve() {

		}
		virtual void Exterior_Derivative_W2V(FaceFieldDv<real, d>& F, const FieldDv<Vector2C, d>& C) {
			Assert(!C.Empty(), "Exterior_Derivative_W2V C->F error: C is empty");
			const real h_bar_over_dx = h_bar / C.grid.dx;
			F.Init(C.grid, MathFunc::Zero<real>());
			const Vector2C* cell = C.Data_Ptr();
			if constexpr (d == 2) C.grid.Exec_Kernel(&W2V_Mapping_Kernel2_Padding0, C.grid, F.Data_Ptr(0), F.Data_Ptr(1), cell, h_bar_over_dx);
			//else if constexpr (d == 3) C.grid.Exec_Kernel(&W2V_Mapping_Kernel3_Padding0, C.grid, F.Data_Ptr(0), F.Data_Ptr(1), F.Data_Ptr(2), cell, h_bar_over_dx);
		}
		virtual void Wave_Function_Correction() {
			wave_function.grid.Exec_Kernel(
				&Wave_Function_Correction_Kernel<d>,
				wave_function.grid,
				wave_function.Data_Ptr(),
				pressure.Data_Ptr(),
				wave_function.grid.dx / h_bar
			);
		}
		virtual void Advance(DriverMetaData& metadata) {
			real dt = metadata.dt;

			//TODO: solve
			Schrödinger_Solve();

			Info("Vel before W2V {}", velocity);
			//wave function to velocity
			Exterior_Derivative_W2V(velocity, wave_function);
			Info("Vel after W2V {}", velocity);

			//projection
			//vel_div=div(velocity)
			ExteriorDerivativePadding0::Apply(vel_div, velocity);
			Info("Div vel {}", vel_div);

			auto [iter, res] = MGPCG.Solve(pressure.Data(), vel_div.Data());
			Info("Solve poisson with {} iters and residual {}", iter, res);
			Info("Pressure {}", pressure);

			//velocity+=grad(p)
			ExteriorDerivativePadding0::Apply(temp_velocity, pressure);
			temp_velocity *= poisson.vol;

			velocity += temp_velocity;
			Info("Vel after projection {}", velocity);

			psi_D.Apply(wave_function);

			ExteriorDerivativePadding0::Apply(vel_div, velocity);
			Info("After projection max div {}", GridEulerFunc::Linf_Norm(vel_div));
			
			Wave_Function_Correction();
		}
	};
}