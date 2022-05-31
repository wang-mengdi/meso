//////////////////////////////////////////////////////////////////////////
// Topology Optimization for Thin Shell
// Copyright (c) (2021-), Fan Feng
// This file is part of CompleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "DiscreteShellTopoOpt.h"
#include "Timer.h"
#include "Common.h"
#include "Optimizer.h"
#include "AuxFunc.h"
#include "ImplicitGeometry.h"
#include <mma/MMASolver.h>

using namespace Meso;
template<int d> 
class DiscreteShellTopologyOptimizer : public Optimizer
{Typedef_VectorD(d); Typedef_MatrixD(d);
public:
	DiscreteShellTopoOpt<d> thin_shell;
	MMASolver* mma_solver;

	//general optimizer variables
	real			obj;								//objective
	Array<real>		grad;								//gradients
	Array<real>		var;								//variables
	Array<real>		var_low_bounds;						//lower bound
	Array<real>		var_up_bounds;						//upper bound
	Array<real>     intmed_var;							//intermediate varaibles
	real			constraint;							//constraint value
	Array<real>		constraint_grads;					//constraint gradients
	
	real target_rho = (real)0.2;						////target average thickness fraction on each vertex
	real avg_rho = (real)0.;							////current average thickness fraction
	real rho_min = (real)1e-3;							
	real rho_max = (real)1;								
	real objective;										////current obj
	int power = 3;										////power index when calculating objective
	real mov_lim = (real)0.1;							////move limit for each iteration, can be tuned 

	////parameters
	real constraint_coef = 1;							////coef for constraints, need to be tuned
	
	Array<real> dobj_drho;									////D_obj/D_rho, rho is the density at each vertex
	Array<real> w;											////inverse hessian times dobj_dx
	Array<real> det_drho;									////dobj_drho without the chain rule part
	Array<real> b_w;										////RHS for calculating w
	SparseMatrix<real> dgrad_drho;							////partial derivative of rho wrt. the gradient of bending+stretching energy
	std::shared_ptr<ImplicitGeometry<d>> target_shape;		////target shape for deformation

	//====================Core optimizer functions=================================
	void Optimize(OptimizerDriverMetaData& meta_data) {
		Compute_Objective();
		Compute_Gradient();
		Compute_Constraint();
		Compute_Constraint_Grad();
		Calculate_Bounds();
		intmed_var = var;				//record the intermediate variables before update

		Assert(std::is_same<real, double>::value, "only double data type is supported by the mma solver");
		mma_solver->Update(ArrayFunc::Data<real, DataHolder::HOST>(var), ArrayFunc::Data<real, DataHolder::HOST>(grad), &constraint, ArrayFunc::Data<real, DataHolder::HOST>(constraint_grads), ArrayFunc::Data<real, DataHolder::HOST>(var_low_bounds), ArrayFunc::Data<real, DataHolder::HOST>(var_up_bounds));
	}

	void Output(const bf::path base_path, const int iter) {
		thin_shell.Output(base_path,iter);
	}

	bool Is_Converged(OptimizerDriverMetaData& meta_data) {
		Array<real> temp = var;
		ArrayFunc::Minus(temp, intmed_var);
		real change = ArrayFunc::Max_Abs<real>(temp);
		
		if (change < meta_data.tol) { return true; }
		return false;
	}

	void Calculate_Bounds() { //may be changed later
		for (int i = 0; i < var.size(); i++) {
			var_up_bounds[i] = std::min(rho_max, var[i] + mov_lim);
			var_low_bounds[i] = std::max(rho_min, var[i] - mov_lim);
		}
	}

	void Init() {
		int vtx_num = thin_shell.Vtx_Num();
		dobj_drho.resize(vtx_num);
		w.resize(vtx_num * d);
		b_w.resize(vtx_num * d);
		det_drho.resize(vtx_num);
		dgrad_drho.resize(vtx_num*d, vtx_num*d);
		thin_shell.Allocate_DGrad_DRho(dgrad_drho);

		var.resize(vtx_num);
		intmed_var.resize(vtx_num);
		grad.resize(vtx_num);
		var_low_bounds.resize(vtx_num);
		var_up_bounds.resize(vtx_num);
		constraint_grads.resize(vtx_num);
		Sync_Var_Fem_To_Opt();
		mma_solver = new MMASolver(vtx_num,1);
	}

	//release memory of the mma solver
	~DiscreteShellTopologyOptimizer() { delete mma_solver; }

	//================== Essential MMA Functions================================

	////rho -> fem_variable_coef -> fem->obj
	virtual real Compute_Objective() {
		thin_shell.Reset_To_Rest_Position();
		Sync_Var_Opt_To_Fem();
		thin_shell.Advance_Quasi_Static(this->H_To_Rho, power);
		Obj();
		return objective;
	}

	////var -> thin_shell.densities -> dobj_drho
	virtual void Compute_Gradient() {
		Update_Grad();
		Sync_Grad_Fem_To_Opt();
		//Numerical_Dobj_Drho();
	}

	////volume constraints
	virtual void Compute_Constraint() {
		// set the total volume fraction is less than frac
		avg_rho = 0;
		for (int i = 0; i < var.size(); i++) { avg_rho += var[i]; }
		avg_rho /= (real)var.size();
		constraint = avg_rho - target_rho;
	}

	virtual void Compute_Constraint_Grad() {
		for (int i = 0; i < constraint_grads.size(); i++) {
			constraint_grads[i] = constraint_coef * (real)1.0 / (real)var.size();
		}
	}

	//====================================================================================

	static inline real Rho_To_H(real rho) {
		return (real) 0.1 * rho; //h=0.1rho
	}

	static inline real Dh_Drho() {
		return (real)0.1;
	}

	static inline real H_To_Rho(real h) {
		return h *(real)10;
	}

	void Sync_Var_Opt_To_Fem() {
#pragma omp parallel for
		for (int i = 0; i < var.size(); i++) { thin_shell.hs[i] = Rho_To_H(var[i]); }
		thin_shell.Update_Lambdas();
	}

	void Sync_Var_Fem_To_Opt() {
#pragma omp parallel for
		for (int i = 0; i < var.size(); i++) { var[i] = H_To_Rho(thin_shell.hs[i]); }
	}

	////rho_grad -> grad
	void Sync_Grad_Fem_To_Opt() {
#pragma omp parallel for
		for (int i = 0; i < var.size(); i++) { grad[i] = dobj_drho[i]; }
	}

	real Obj() {
		objective = 0;
		if (target_shape) {
			objective += Target_Deformation_Obj();
		}
		else {
			objective += /*thin_shell.Total_Bending_Obj(this->H_To_Rho,power) +*/ thin_shell.Total_Stretching_Obj(this->H_To_Rho,power);
		}
		return objective;
	}

	real Target_Deformation_Obj() {
		real deformation_obj = 0;
		if (target_shape) {
			for (int i = 0; i < thin_shell.Vtx_Num(); i++) {
				real phi = target_shape->Phi(thin_shell.X()[i]);
				deformation_obj += abs(phi);
			}
		}
		return deformation_obj;
	}

	void Update_DDeformation_DU(Array<real>& dEt_dU) {
		if (target_shape) {
			for (int i = 0; i < thin_shell.Vtx_Num(); i++) {
				real phi = target_shape->Phi(thin_shell.X()[i]);
				VectorD normal = target_shape->Normal(thin_shell.X()[i]);
				if (phi < 0) { normal = -normal; }
				for (int axis = 0; axis < d; axis++) {
					dEt_dU[i * d + axis] +=  normal[axis];
				}
			}
		}
	}

	////gradients
	virtual void Update_Grad() {
		ArrayFunc::Fill(det_drho, 0);
		ArrayFunc::Fill(b_w, 0);

		if (!target_shape) {
			thin_shell.Update_DEt_Stretch_DRho(det_drho, this->Dh_Drho,this->H_To_Rho, power);
			/*thin_shell.Update_DEt_Bend_DRho(det_drho, this->Dh_Drho, this->H_To_Rho, power);*/
			thin_shell.Update_W(b_w, w, dgrad_drho, true, this->Dh_Drho, this->H_To_Rho, power);
		}
		else {
			//Solve for W and update dgrad_drho at the same time
			//if only want the deformation target
			Update_DDeformation_DU(b_w);
			thin_shell.Update_W(b_w, w, dgrad_drho, false, this->Dh_Drho, this->H_To_Rho, power);
		}

		SparseMatrixMapping<real, DataHolder::DEVICE> meso_mat(dgrad_drho); //This maybe changed later. Do it on host?
		ArrayDv<real> Dv_tmp;
		ArrayDv<real> Dv_w = w;
		meso_mat.Apply(Dv_tmp, Dv_w);
		Array<real>tmp = Dv_tmp;
		//Get dL/drho
#pragma omp parallel for
		for (int i = 0; i < thin_shell.Vtx_Num(); i++) {dobj_drho[i] = det_drho[i] - tmp[i*d];} //extract the only useful part from the padded tmp
	}


	////=======================numerical verification, not sure how to do it correctly=========================================
	void Numerical_Dobj_Drho() {
		real old_obj = objective;
		real drho = (real)1e-7;
		Array<real> numeric_dobj_drho(dobj_drho.size(), (real)0);

		std::cout << "Numerical derivative for dobj_drho on each vertex" << std::endl;
		for (int i = 0; i < var.size(); i++) {
			real old_h = thin_shell.hs[i];
			thin_shell.Reset_To_Rest_Position();
			thin_shell.hs[i] += Dh_Drho()*drho;
			thin_shell.Update_Lambdas();
			thin_shell.Advance_Quasi_Static(this->H_To_Rho,power);
			real obj2 = Obj();
			numeric_dobj_drho[i] = (obj2 - old_obj) / drho;
			if (abs(numeric_dobj_drho[i] - dobj_drho[i]) / abs(numeric_dobj_drho[i]) > (real)5e-2) { // print when the difference is big
				std::cout << "vertex " << i << ": rho ["<< H_To_Rho(thin_shell.hs[i]) <<"], ana [" << dobj_drho[i] << "], num [" << numeric_dobj_drho[i] << "]" << std::endl;
			}
			thin_shell.hs[i] = old_h;
		}
	}
};

