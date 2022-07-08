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
#include "Advection.h"

namespace Meso {

	template<int d>
	class FluidPIC : public Simulator {
		Typedef_VectorD(d);
	public:
		FieldDv<real, d> pressure;
		FieldDv<real, d> vel_div;
		FaceField<real, d> weight;
		FaceField<real, d> velocity_host;
		FaceFieldDv<real, d> velocity_dev;
		FaceFieldDv<real, d> temp_velocity;

		BoundaryConditionDirect<FaceFieldDv<real, d>> psi_N;
		MaskedPoissonMapping<real, d> poisson;
		VCycleMultigridIntp<real, d> MG_precond;
		ConjugateGradient<real> MGPCG;
		
		Array<VectorD> particles_pos;
		Array<VectorD> particles_vel;

		void Clamp_Particles_Pos() {
			Grid<d> grid = velocity_host.grid;
			VectorD pos_min = grid.Node_Min();
			VectorD pos_max = grid.Node_Max();
			const int particles_num = particles_pos.size();
			for (int p = 0; p < particles_num; p++)
				for (int axis = 0; axis < d; axis++)
					particles_pos[p][axis]= std::min(pos_max(axis), std::max(pos_min(axis), particles_pos[p][axis]));
		}

		void Reinit_Particles(){
			particles_pos.clear();
			particles_vel.clear();
			Grid<d> grid = velocity_host.grid;
			grid.Iterate_Nodes([&](const VectorDi& cell) {
				const VectorD cell_center_pos = grid.Position(cell);
				const int particles_per_cell = (1 << d);
				for (int i = 0; i < particles_per_cell; i++)
				{
					VectorD particle_pos = cell_center_pos + VectorD::Random() * grid.dx * 0.5;
					particles_pos.push_back(particle_pos);
					particles_vel.push_back(VectorD::Zero());
				}
				});
			Clamp_Particles_Pos();
		}

		void Init(Field<bool, d>& fixed, FaceField<real, d>& vol, FaceField<bool, d>& face_fixed, FaceField<real, d>& initial_velocity) {
			Grid<d> grid = fixed.grid;
			pressure.Init(grid);
			vel_div.Init(grid);
			weight.Init(grid);
			velocity_host.Deep_Copy(initial_velocity);
			velocity_dev.Deep_Copy(initial_velocity);
			temp_velocity.Deep_Copy(initial_velocity);
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

		void Grid_To_Particles(){
			const int particles_num = particles_pos.size();
			for (int p = 0; p < particles_num; p++)
				particles_vel[p] = IntpLinear::Face_Vector(velocity_host, particles_pos[p]);
			Clamp_Particles_Pos();
		}

		void Particles_Operation(const real dt){
			const int particles_num = particles_pos.size();
			for (int p = 0; p < particles_num; p++)
				particles_pos[p] += particles_vel[p] * dt;
		}

		void Particles_To_Grid(){
			weight.Fill(0);
			velocity_host.Fill(0);
			
			static constexpr int dx[8] = { 0,1,0,1,0,1,0,1 };
			static constexpr int dy[8] = { 0,0,1,1,0,0,1,1 };
			static constexpr int dz[8] = { 0,0,0,0,1,1,1,1 };
			const int particles_num = particles_pos.size();
			for (int p = 0; p < particles_num; p++) 
				for (int axis = 0; axis < d; axis++) {
					Grid<d> face_grid = velocity_host.grid.Face_Grid(axis);
					VectorDi node; VectorD frac;
					face_grid.Get_Fraction(particles_pos[p], node, frac);
					VectorDi face_counts = face_grid.Counts();
					if constexpr (d == 2) {
						real w[2][2] = { {1.0 - frac[0],frac[0]},{1.0 - frac[1],frac[1]} };
						for (int s = 0; s < 4; s++)
						{
							int d0 = dx[s], d1 = dy[s];
							int idx = face_grid.Index(node[0] + d0, node[1] + d1);
							(*velocity_host.face_data[axis])[idx] += w[0][d0] * w[1][d1] * particles_vel[p][axis];
							(*weight.face_data[axis])[idx] += w[0][d0] * w[1][d1];
						}
					}
					if constexpr (d == 3) {
						real w[3][2] = { {1.0 - frac[0],frac[0]},{1.0 - frac[1],frac[1]} ,{1.0 - frac[2],frac[2]} };
						for (int s = 0; s < 8; s++)
						{
							int d0 = dx[s], d1 = dy[s], d2 = dz[s];
							int idx = face_grid.Index(node[0] + d0, node[1] + d1, node[2] + d2);
							(*velocity_host.face_data[axis])[idx] += w[0][d0] * w[1][d1] * w[2][d2] * particles_vel[p][axis];
							(*weight.face_data[axis])[idx] += w[0][d0] * w[1][d1] * w[2][d2];
						}
					}
				}
			velocity_host.Calc_Faces([&](const int axis, const VectorDi& face) ->real{
				int idx = weight.grid.Face_Grid(axis).Index(face);
				if (( * weight.face_data[axis])[idx] > 0)
					return (* velocity_host.face_data[axis])[idx] / (*weight.face_data[axis])[idx];
				else
					return 0;
				});
			velocity_dev.Deep_Copy(velocity_host);
		}

		void Grid_Operation(){
			psi_N.Apply(velocity_dev);

			ExteriorDerivativePadding0::Apply(vel_div, velocity_dev);

			auto [iter, res] = MGPCG.Solve(pressure.Data(), vel_div.Data());
			Info("Solve poisson with {} iters and residual {}", iter, res);

			ExteriorDerivativePadding0::Apply(temp_velocity, pressure);
			temp_velocity *= poisson.vol;

			velocity_dev += temp_velocity;
			psi_N.Apply(velocity_dev);

			ExteriorDerivativePadding0::Apply(vel_div, velocity_dev);
			Info("After projection max div {}", GridEulerFunc::Linf_Norm(vel_div));
			
			velocity_host.Deep_Copy(velocity_dev);
		}

		virtual void Advance(DriverMetaData& metadata) {
			Grid_To_Particles();
			Particles_Operation(metadata.dt);
			Particles_To_Grid();
			Grid_Operation();
			Reinit_Particles();
		}
	};
}