//////////////////////////////////////////////////////////////////////////
// Fluid Euler
// Copyright (c) (2022-), Bo Zhu, Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once

#include "Multigrid.h"
#include "ConjugateGradient.h"
#include "Advection.h"
#include "Simulator.h"
#include "IOFunc.h"
#include "GridEulerFunc.h"

namespace Meso {
	template<int d>
	class FluidEuler : public Simulator {
		Typedef_VectorD(d);
	public:
		FaceFieldDv<real, d> velocity;
		BoundaryConditionDirect<FaceFieldDv<real, d>> psi_N;

		FaceFieldDv<real, d> temp_velocity;
		FieldDv<real, d> pressure;
		FieldDv<real, d> vel_div;
		MaskedPoissonMapping<real, d> poisson;
		VCycleMultigridIntp<real, d> MG_precond;
		ConjugateGradient<real, d> MGPCG;

		FaceField< real, d> velocity_host;
		Field<real, d> vorticity;

		void Init(Field<unsigned char, d>& cell_type, FaceField<real, d>& vol, FaceField<bool, d>& face_fixed, 
			FaceField<real, d>& initial_velocity, bool is_pure_neumann = false) {
			velocity.Deep_Copy(initial_velocity);
			psi_N.Init(face_fixed, initial_velocity);
			temp_velocity.Init(velocity.grid);
			pressure.Init(velocity.grid);
			vel_div.Init(velocity.grid);
	
			poisson.Init(cell_type, vol);
			MG_precond.Init_Poisson(poisson, 2);
			MGPCG.Init(&poisson, &MG_precond, false, -1, 1e-5, is_pure_neumann);

			Grid<d> vel_grid = velocity.grid;
			velocity_host.Init(vel_grid);
			Grid<d> vor_grid(vel_grid.Counts() + VectorDi::Ones(), vel_grid.dx, vel_grid.Domain_Min(CENTER), CORNER);
			vorticity.Init(vor_grid);
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

			Output_Vorticity(metadata);
		}
		virtual void Advance(DriverMetaData& metadata) {
			real dt = metadata.dt;

			//advection
			Advection<IntpLinearClamp, 2>::Advect(dt, temp_velocity, velocity, velocity);
			velocity = temp_velocity;
			psi_N.Apply(velocity);

			//projection
			//vel_div=div(velocity)
			ExteriorDerivativePadding0::Apply(vel_div, velocity);

			auto [iter, res] = MGPCG.Solve(pressure.Data(), vel_div.Data());
			Info("Solve poisson with {} iters and residual {}", iter, res);

			//velocity+=grad(p)
			ExteriorDerivativePadding0::Apply(temp_velocity, pressure);
			temp_velocity *= poisson.vol;
			velocity += temp_velocity;

			ExteriorDerivativePadding0::Apply(vel_div, velocity);
			Info("After projection max div {}", GridEulerFunc::Linf_Norm(vel_div));
		}

		virtual void Output_Vorticity(DriverMetaData& metadata) {
			if constexpr (d == 2) {
				//calculate vorticity
				velocity_host = velocity;
				vorticity.Calc_Nodes([&] (const VectorDi cell) {
					for (int axis = 0; axis < 2; axis++)
						if (cell[axis] == 0 || cell[axis] == vorticity.grid.Counts()[axis] - 1)
							return (real)0;
				return -(velocity_host(0, cell) - velocity_host(0, cell - VectorDi::Unit(1)) + 
					velocity_host(1, cell - VectorDi::Unit(0)) - velocity_host(1, cell)) / vorticity.grid.dx;
				});
				//output vorticity
				std::string vorticity_name = fmt::format("vor{:04d}.vts", metadata.current_frame);
				bf::path vorticity_path = metadata.base_path / bf::path(vorticity_name);
				VTKFunc::Write_VTS(vorticity, vorticity_path.string());
			}
		}
	};
}