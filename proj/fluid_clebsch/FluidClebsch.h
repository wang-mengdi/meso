//////////////////////////////////////////////////////////////////////////
// Clebsch Fluid: Implementation of Schrödinger's Smoke [Chern et al. 2015]
// Copyright (c) (2022-), Bo Zhu, Mengdi Wang, Zhecheng Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once

#include <array>

#include "Multigrid.h"
#include "ConjugateGradient.h"
#include "ExteriorDerivative.h"
#include "Advection.h"
#include "Simulator.h"
#include "IOFunc.h"
#include "AuxFunc.h"
#include "GridEulerFunc.h"
#include "WaveFunc.h"
#include "device_launch_parameters.h"
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>

namespace Meso {
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
		ConjugateGradient<real> CG;
		void Init(real h_bar_, Field<bool, d>& fixed, Field<Vector2C, d>& initial_wave_function, FaceField<real, d>& vol) {
			h_bar = h_bar_;
			wave_function.Deep_Copy(initial_wave_function);
			WaveFunc::Exterior_Derivative_W2V(velocity, wave_function, h_bar);
			psi_D.Init(fixed, initial_wave_function);
			temp_velocity.Init(velocity.grid);
			pressure.Init(velocity.grid);
			vel_div.Init(velocity.grid);
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
		virtual void Schrodinger_Solve() {

		}
		virtual void Normalize() {
			wave_function.grid.Exec_Kernel(
				&WaveFunc::Wave_Function_Normalization_Kernel<d>,
				wave_function.grid,
				wave_function.Data_Ptr()
			);
		}

		virtual void Wave_Function_Correction() {
			wave_function.grid.Exec_Kernel(
				&WaveFunc::Wave_Function_Correction_Kernel<d>,
				wave_function.grid,
				wave_function.Data_Ptr(),
				pressure.Data_Ptr(),
				wave_function.grid.dx / h_bar
			);
		}
		virtual void Advance(DriverMetaData& metadata) {
			real dt = metadata.dt;

			//TODO: solve 
			Schrodinger_Solve();

			Normalize();

			Info("Vel before W2V {}", velocity);
			//wave function to velocity
			WaveFunc::Exterior_Derivative_W2V(velocity, wave_function, h_bar);
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

			ExteriorDerivativePadding0::Apply(vel_div, velocity);
			Info("After projection max div {}", GridEulerFunc::Linf_Norm(vel_div));
			
			Wave_Function_Correction();
			psi_D.Apply(wave_function);
		}
	};
}