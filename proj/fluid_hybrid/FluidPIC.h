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

		int Delete_Particles(const Array<bool>& is_deleted){
			const int particles_num = particles_pos.size();
			int i = 0;
			for (int p = 0; p < particles_num; p++)
				if (!is_deleted[p])
				{
					particles_pos[i] = particles_pos[p];
					particles_vel[i] = particles_vel[p];
					++i;
				}
			particles_pos.resize(i);
			particles_vel.resize(i);
			return i;
		}

		bool Out_Of_Range(const VectorD& particle_pos, const VectorD&pos_min, const VectorD& pos_max) {
			for (int axis = 0; axis < d; axis++)
				if (particle_pos[axis]<pos_min[axis] || particle_pos[axis]>pos_max[axis])
					return true;
			return false;
		}

		void Delete_Out_Of_Range_Particles(){
			const int particles_num = particles_pos.size();
			Array<bool> is_deleted(particles_num, false);
			VectorD pos_min = velocity_host.grid.Node_Min();
			VectorD pos_max = velocity_host.grid.Node_Max();
			for (int p = 0; p < particles_num; p++)
				is_deleted[p] = Out_Of_Range(particles_pos[p], pos_min, pos_max);
			Delete_Particles(is_deleted);
		}

		void Init_Particles(){

			Grid<d> grid = velocity_host.grid;
			const int particles_per_cell = (1 << d);
			const int particles_num = grid.Counts().prod() * particles_per_cell;
			particles_pos.resize(particles_num);
			particles_vel.resize(particles_num);

			int p = 0;
			VectorD pos_min = grid.Node_Min();
			VectorD pos_max = grid.Node_Max();
			grid.Iterate_Nodes([&](const VectorDi& cell) {
				const VectorD cell_center_pos = grid.Position(cell);
				for (int i = 0; i < particles_per_cell; i++)
				{
					VectorD delta_pos= VectorD::Random() * grid.dx * 0.5;
					while (Out_Of_Range(cell_center_pos + delta_pos, pos_min, pos_max))
						delta_pos = VectorD::Random() * grid.dx * 0.5;
					particles_pos[p] = cell_center_pos + delta_pos;
					particles_vel[p] = VectorD::Zero();
					++p;
				}
				});
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
			Grid_Operation();
			
			Init_Particles();
		}

		virtual real CFL_Time(const real cfl) {
			real dx = velocity_host.grid.dx;
			real max_vel = GridEulerFunc::Linf_Norm(velocity_host);
			return dx * cfl / max_vel;
		}

		virtual void Output_Vorticity(DriverMetaData& metadata) {
			if constexpr (d == 2) {
				//calculate vorticity
				Field<real, d> vorticity;
				Grid<d> grid = velocity_host.grid;
				vorticity.Init(grid);
				vorticity.Calc_Nodes([&](const VectorDi cell) {
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

		virtual void Output(DriverMetaData& metadata) {
			std::string grid_name = fmt::format("grid{:04d}.vts", metadata.current_frame);
			bf::path grid_path = metadata.base_path / bf::path(grid_name);
			VTKFunc::Write_VTS(velocity_host, grid_path.string());
			Output_Vorticity(metadata);
			std::string particles_name = fmt::format("particles{:04d}.vtu", metadata.current_frame);
			bf::path particles_path = metadata.base_path / bf::path(particles_name);
			VTKFunc::Write_VTU_Particles<d>(particles_pos, particles_vel, particles_path.string());
		}

		void Reseeding() {
			Grid<d> grid = velocity_host.grid;
			Field<int, d> particles_num_in_cell;
			particles_num_in_cell.Init(grid, 0);
			int particles_num = particles_pos.size();
			Array<VectorDi> particle_to_cell(particles_num);
			
			VectorD pos_min = grid.Node_Min();
			real dx = grid.dx;
			for (int p = 0; p < particles_num; p++)
			{
				VectorDi cell = ((particles_pos[p] - pos_min) / dx).template cast<int>();
				particle_to_cell[p] = cell;
				++particles_num_in_cell(cell);
			}
			
			Array<bool> is_deleted(particles_num, false);
			const int count_max = (1 << (d + 1));
			for (int p = 0; p < particles_num; p++)
				if (particles_num_in_cell(particle_to_cell[p]) > count_max){
					is_deleted[p] = true;
					--particles_num_in_cell(particle_to_cell[p]);
				}
			particles_num = Delete_Particles(is_deleted);

			int add_particles_num = 0;
			particles_num_in_cell.Iterate_Nodes([&](const VectorDi cell) {
				if (particles_num_in_cell(cell) < 2)
					add_particles_num += 2 - particles_num_in_cell(cell);
				});
			
			particles_pos.resize(particles_num + add_particles_num);
			particles_vel.resize(particles_num + add_particles_num);
			VectorD pos_max = grid.Node_Max();
			particles_num_in_cell.Iterate_Nodes([&](const VectorDi cell) {
				while (particles_num_in_cell(cell) < 2){
					VectorD cell_center_pos = grid.Position(cell);
					VectorD delta_pos = VectorD::Random() * grid.dx * 0.5;
					while (Out_Of_Range(cell_center_pos + delta_pos, pos_min, pos_max))
						delta_pos = VectorD::Random() * grid.dx * 0.5;
					particles_pos[particles_num] = cell_center_pos + delta_pos;
					particles_vel[particles_num] = VectorD::Zero();
					++particles_num;
					++particles_num_in_cell(cell);
				}
				});

		}

		void Grid_To_Particles(){
			const int particles_num = particles_pos.size();
			for (int p = 0; p < particles_num; p++)
				particles_vel[p] = IntpLinearClamp::Face_Vector(velocity_host, particles_pos[p]);
		}

		void Particles_Operation(const real dt){
			const int particles_num = particles_pos.size();
			for (int p = 0; p < particles_num; p++)
				particles_pos[p] += particles_vel[p] * dt;
			Delete_Out_Of_Range_Particles();
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
					//clamp
					for (int axis = 0; axis < d; axis++) {
						if (node[axis] < 0) node[axis] = 0, frac[axis] = 0;
						if (node[axis] > face_counts[axis] - 2) node[axis] = face_counts[axis] - 2, frac[axis] = 1;
					}
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
			Reseeding();
			Grid_To_Particles();
			Particles_Operation(metadata.dt);
			Particles_To_Grid();
			Grid_Operation();
		}
	};
}