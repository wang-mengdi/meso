#pragma once
#include "FluidEuler.h"
#include "Multigrid.h"
#include "ConjugateGradient.h"
#include "Advection.h"
#include "Stretching.h"
#include "Simulator.h"
#include "IOFunc.h"
namespace Meso {
	template<int d>
	class FluidImpulse : public Simulator {
		Typedef_VectorD(d);
	public:
		FaceFieldDv<real, d> velocity;
		FieldDv<VectorD, d> inverse_flow_map;
		FaceFieldDv<VectorD, d> inverse_flow_map_grad;
		BoundaryConditionDirect<FaceFieldDv<real, d>> psi_N;

		FaceFieldDv<real, d> temp_velocity;
		FieldDv<real, d> pressure;
		FieldDv<real, d> vel_div;
		MaskedPoissonMapping<real, d> poisson;
		VCycleMultigridIntp<real, d> MG_precond;
		ConjugateGradient<real> MGPCG;

		void Init(Field<bool, d>& fixed, FaceField<real, d>& vol, FaceField<bool, d>& face_fixed, FaceField<real, d>& initial_velocity) {
			velocity.Deep_Copy(initial_velocity);
			inverse_flow_map.Init(velocity.grid);
			inverse_flow_map_grad.Init(velocity.grid);
			psi_N.Init(face_fixed, initial_velocity);
			temp_velocity.Init(velocity.grid);
			pressure.Init(velocity.grid);
			vel_div.Init(velocity.grid);

			poisson.Init(velocity.grid, vol, fixed);
			MG_precond.Init_Poisson(poisson, 2, 2);
			MGPCG.Init(&poisson, &MG_precond, false, -1, 1e-6);
		}

		virtual real CFL_Time(const real cfl) {
			real dx = velocity.grid.dx;
			real max_vel = velocity.Max_Abs();
			return dx * cfl / max_vel;
		}

		virtual void Output(const bf::path base_path, const int frame) {
			std::string vts_name = fmt::format("vts{:04d}.vts", frame);
			bf::path vtk_path = base_path / bf::path(vts_name);
			VTKFunc::Write_VTS(velocity, vtk_path.string());
		}

		virtual void Advance(const int current_frame, const real current_time, const real dt) {
			//advection
			SemiLagrangian::Advect(dt, temp_velocity, velocity, velocity);
			SemiLagrangian::Inverse_Flow_Map(dt,inverse_flow_map,velocity);
			velocity = temp_velocity;
			
			//stretching
			Exterior_Derivative<VectorD,d>(inverse_flow_map_grad, inverse_flow_map);
			inverse_flow_map_grad *= (real)1 / velocity.grid.dx;
			Stretching::Covector_Stretching(temp_velocity, inverse_flow_map_grad, velocity);
			velocity = temp_velocity;	//Not sure with the boundary condition here
			psi_N.Apply(velocity);
			//projection
			Exterior_Derivative(vel_div, velocity);

			int iter; real res;
			MGPCG.Solve(pressure.Data(), vel_div.Data(), iter, res);
			Info("Solve poisson with {} iters and residual {}", iter, res);

			Exterior_Derivative(temp_velocity, pressure);
			temp_velocity *= poisson.vol;

			velocity += temp_velocity;
			psi_N.Apply(velocity);

			Exterior_Derivative(vel_div, velocity);
			Info("After projection max div {}", vel_div.Max_Abs());
		}
	};
}