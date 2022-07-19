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
		ConjugateGradient<real> MGPCG;
		void Init(Field<bool, d>& fixed, FaceField<real, d>& vol, FaceField<bool, d>& face_fixed, FaceField<real, d>& initial_velocity) {
			velocity.Deep_Copy(initial_velocity);
			psi_N.Init(face_fixed, initial_velocity);
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

			Output_Vorticity(metadata);
		}
		virtual void Advance(DriverMetaData& metadata) {
			real dt = metadata.dt;

			//advection
			SemiLagrangian<IntpLinearPadding0>::Advect(dt, temp_velocity, velocity, velocity);
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
			psi_N.Apply(velocity);

			ExteriorDerivativePadding0::Apply(vel_div, velocity);
			Info("After projection max div {}", GridEulerFunc::Linf_Norm(vel_div));
		}

		virtual void Output_Vorticity(DriverMetaData& metadata) {
			if constexpr (d == 2) {
				//calculate vorticity
				FaceField < real, d> velocity_host = velocity;
				Field<real, d> vorticity;
				Grid<d> grid = velocity_host.grid;
				vorticity.Init(grid);
				vorticity.Calc_Nodes([&] (const VectorDi cell) {
					for (int axis = 0; axis < 2; axis++)
						if (cell[axis] == 0 || cell[axis] == grid.Counts()[axis] - 1)
							return (real)0;
					real val[2];
					for (int axis = 0; axis < 2; axis++) {
						VectorDi face_0 = cell + VectorDi::Ones();
						VectorDi face_1 = cell + VectorDi::Unit(axis);
						VectorDi face_2 = cell - VectorDi::Unit(axis);
						VectorDi face_3 = cell - VectorDi::Unit(axis) + VectorDi::Unit(1 - axis);
						val[axis] = velocity_host.Get(1 - axis, face_0) + velocity_host.Get(1 - axis, face_1);
						val[axis] -= velocity_host.Get(1 - axis, face_2) + velocity_host.Get(1 - axis, face_3);
						val[axis] *= 0.25 / grid.dx;
					}
					return val[0] - val[1];
				});
				//output vorticity
				std::string vorticity_name = fmt::format("vor{:04d}.vts", metadata.current_frame);
				bf::path vorticity_path = metadata.base_path / bf::path(vorticity_name);
				VTKFunc::Write_VTS(vorticity, vorticity_path.string());
			}
		}
	};
}