//////////////////////////////////////////////////////////////////////////
// Fluid PIC
// Copyright (c) (2022-), Yuchen Sun
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once

#include "Simulator.h"
#include "Multigrid.h"
#include "ConjugateGradient.h"
#include "GridEulerFunc.h"
#include "IOFunc.h"

namespace Meso {

	template<int d>
	class FluidPIC : public Simulator {
		Typedef_VectorD(d);
	public:
		FieldDv<real, d> pressure;
		FieldDv<real, d> vel_div;
		FaceField<real, d> velocity_host;
		FaceField<real, d> weight_host;
		FaceFieldDv<real, d> velocity_dev;
		FaceFieldDv<real, d> weight_dev;

		BoundaryConditionDirect<FaceFieldDv<real, d>> psi_N;
		MaskedPoissonMapping<real, d> poisson;
		VCycleMultigridIntp<real, d> MG_precond;
		ConjugateGradient<real> MGPCG;
		
		Array<VectorD> particles_pos;
		Array<VectorD> particles_vel;

		void Reinit_Particles()
		{
			particles_pos.clear();
			particles_vel.clear();
			Grid<d> grid = velocity_host.grid;
			grid.Iterate_Nodes([&](const VectorDi& cell) {
				const VectorD cell_center_pos = grid.Position(cell);
				for (int i = 0; i < (1 << d); i++)
				{
					VectorD particle_pos = cell_center_pos + VectorD::Random() * grid.dx * 0.5;
					particles_pos.push_back(particle_pos);
					particles_vel.push_back(VectorD::Ones());
				}
				});
		}

		void Init(Field<bool, d>& fixed, FaceField<real, d>& vol, FaceField<bool, d>& face_fixed, FaceField<real, d>& initial_velocity) {
			Grid<d> grid = fixed.grid;
			pressure.Init(grid);
			vel_div.Init(grid);
			velocity_host.Deep_Copy(initial_velocity);
			weight_host.Init(fixed.grid, 0);
			velocity_dev.Deep_Copy(initial_velocity);
			weight_dev.Init(fixed.grid, 0);
			psi_N.Init(face_fixed, initial_velocity);
			poisson.Init(fixed, vol);
			MG_precond.Init_Poisson(poisson, 2, 2);
			MGPCG.Init(&poisson, &MG_precond, false, -1, 1e-6);
			Reinit_Particles();
		}

		virtual real CFL_Time(const real cfl) {
			real dx = velocity_host.grid.dx;
			real max_vel = GridEulerFunc::Linf_Norm(velocity_host);
			return dx * cfl / max_vel;
		}

		virtual void Output(DriverMetaData& metadata) {
			std::string grid_name = fmt::format("grid{:04d}.vts", metadata.current_frame);
			bf::path grid_path = metadata.base_path / bf::path(grid_name);
			VTKFunc::Write_VTS(velocity_host, grid_path.string());
			std::string particles_name = fmt::format("particles{:04d}.vtu", metadata.current_frame);
			bf::path particles_path = metadata.base_path / bf::path(particles_name);
			VTKFunc::Write_VTU_Particles<d>(particles_pos, particles_vel, particles_path.string());
		}

		virtual void Advance(DriverMetaData& metadata) {
		}
	};
}