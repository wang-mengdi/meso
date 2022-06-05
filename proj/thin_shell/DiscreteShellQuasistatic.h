#pragma once
#include "DiscreteShell.h"
#include "Optimizer.h"
#include "IOHelper.h"
#include "AuxFunc.h"

using namespace Meso;
enum class OptimizeMode { NewtonRaphson, LevenbergMarquardt };

template<int d> class DiscreteShellQuasistatic:public DiscreteShell<d>, public Optimizer
{
	Typedef_VectorD(d); Typedef_MatrixD(d); using Base = DiscreteShell<d>;

public:
	real err = (real)1; real initial_potential_energy;
	OptimizeMode optimize_mode=OptimizeMode::LevenbergMarquardt;

	void Output(const bf::path base_path, const int iter) {
		std::string vtu_name = fmt::format("vtu{:04d}.vtu", iter);
		bf::path vtu_path = base_path / bf::path(vtu_name);
		DiscreteShellVTKFunc::Output_VTU<d, VectorD>(mesh, V(), vtu_path.string());
	}

	void Optimize(OptimizerDriverMetaData& meta_data) //optimize for displacement dx to minimize energy
	{
		switch (optimize_mode) {
		case OptimizeMode::NewtonRaphson:
			Optimize_NewtonRaphson(meta_data);
			break;
		case OptimizeMode::LevenbergMarquardt:
			Optimize_LevenbergMarquardt(meta_data);
			break;
		default:
			Error("Only NewtonRaphson and LevenbergMarquardt modes are available.");
			break;
		}
	}

	void Optimize_NewtonRaphson(OptimizerDriverMetaData& meta_data) {
		damping = (real)0;

		Base::Clear_A_And_Rhs();

		//add external forces
		for (auto& force : bc.forces) {
			Base::Add_Block(b, force.first, force.second);
		}

		Base::Update_Implicit_Stretching((real)1.0); //equivalent to dt=1.0
		Base::Update_Implicit_Bending((real)1.0);
		Base::Update_Implicit_Boundary_Condition();
		Base::Solve();

		Update_Displacement();

		Info("Total energy: {}", Total_Potential_Energy());
		err = ArrayFunc::Norm<real, DataHolder::HOST>(dx) / dx.size();
	}

	void Optimize_LevenbergMarquardt(OptimizerDriverMetaData& meta_data) {
		if (meta_data.iter_count == 1) {initial_potential_energy = Total_Potential_Energy();}

		real mu = 0; real nu = 2;

		damping = (real)0;

		Base::Clear_A_And_Rhs();

		//add external forces
		for (auto& force : bc.forces) {
			Base::Add_Block(b, force.first, force.second);
		}

		Base::Update_Implicit_Stretching((real)1.0); //equivalent to dt=1.0
		Base::Update_Implicit_Bending((real)1.0);
		Base::Update_Implicit_Boundary_Condition();

		//K=K+muI
		SparseMatrix<real> I(Base::Vtx_Num()*d, Base::Vtx_Num()*d);
		I.setIdentity();
		Base::A+=mu * I;

		Base::Solve();

		Update_Displacement();

		Info("Total energy: {}", Total_Potential_Energy());
		err = ArrayFunc::Norm<real, DataHolder::HOST>(dx) / dx.size();
	}

	void Update_Displacement() {
#pragma omp parallel for
		for (int i = 0; i < Vtx_Num(); i++) {
			for (int j = 0; j < d; j++) {
				X()[i][j] += dx[i * d + j]; //may multiply stepsize
			}
		}
	}

	bool Is_Converged(OptimizerDriverMetaData& meta_data) {
		if (err < meta_data.tol) { return true; }
		else { return false; }
	}

	real External_Force_Potential_Energy() {
		real potential_energy = (real)0;
		for (auto& force : bc.forces) {
			potential_energy -= force.second.dot(X()[force.first]); //potential energy by the external force
		}
		return potential_energy;
	}

	real Total_Potential_Energy() {
		return Base::Total_Stretching_Energy() + Base::Total_Bending_Energy() + External_Force_Potential_Energy();
	}

	real Quadratic_Energy(real totoal_energy) {
		return totoal_energy + (real)1.5*ArrayFunc::Dot(Base::dx, Base::b);
	}
};