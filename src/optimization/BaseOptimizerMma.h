//#pragma once
//
//#include "Common.h"
//#include "AuxFunc.h"
//#include "Timer.h"
//#include "MetaData.h"
//#include <mma/MMASolver.h>
//#include <iostream>
//
//using namespace Meso;
//class OptimizerMMA
//{
//public:
//	int n_var;			////var num
//	int n_cons;			////constraints num
//	bool use_var_L;		////lower bound
//	bool use_var_U;		////upper bound
//
//	////bound values
//	real var_lb;		////lower bound for var
//	real var_ub;		////upper bound for var
//	real cons_lb;		////lower bound for constraints
//	real cons_ub;		////upper bound for constraints
//
//	//used to store intermediate values
//	real intmed_obj;
//	Array<real> intmed_var;		//written depend on write_intmed
//
//	real* var = nullptr;		////variable
//	real* var_L = nullptr;		////current lower bound
//	real* var_U = nullptr;		////current upper bound
//
//	bool use_elementwise_bound = false;
//	real* var_LB = nullptr;		////preset elementwise lower bound
//	real* var_UB = nullptr;		////preset elementwise upper bound
//
//	real* grad = nullptr;
//	real* cons = nullptr;
//	real* cons_grad = nullptr;
//
//	real obj = (real)std::numeric_limits<real>::max();
//	real movlim = (real).2;
//	real* movlims = nullptr;		////preset movlimits, used when using elementwise bound
//
//	OptimizerMMA() {}
//	~OptimizerMMA() { Delete_Data(); }
//
//	//needs to be realized in the child class
//	virtual real Compute_Objective(const real* var) { return (real)0; }
//	virtual void Compute_Gradient(const real* var, real* grad) {}
//	virtual void Compute_Constraint(const real* var, real* constraints) { if (n_cons == 0)return; }
//	virtual void Compute_Constraint_Grad(const real* var, real* constraint_grads) {}
//	virtual void Write_Substep(const int frame) {}
//
//protected:
//	template<class T_VAL> T_VAL* New(const int n) { if (n == 0)return nullptr; else return (T_VAL*)malloc(sizeof(T_VAL) * n); }
//	template<class T_VAL> void Delete(T_VAL*& ptr) { if (ptr != nullptr) { free(ptr); ptr = nullptr; } }
//
//	void Allocate_Data()
//	{
//		////assuming n_cons,n_var already initialized
//		Delete(var); var = New<real>(n_var);
//		Delete(var_L); var_L = New<real>(n_var);
//		Delete(var_U); var_U = New<real>(n_var);
//
//		if (use_elementwise_bound) {
//			Delete(var_LB); var_LB = New<real>(n_var);
//			Delete(var_UB); var_UB = New<real>(n_var);
//			Delete(movlims); movlims = New<real>(n_var);
//		}
//
//		Delete(grad); grad = New<real>(n_var);
//		Delete(cons); cons = New<real>(n_cons);
//		Delete(cons_grad); cons_grad = New<real>(n_var * n_cons);
//
//		intmed_var.resize(n_var);
//	}
//
//	virtual void Optimize(json& j)
//	{
//		OptimizerMetaData meta_data;
//		meta_data.Init(j);
//
//		MMASolver* mma = new MMASolver(n_var, n_cons);
//		//mma->ConstraintModification(true); //Fan: what does this mean? Ihe function is never used
//
//		////assuming var_LB and var_UB were set outside the class
//		real ch = 1.0;
//		while (ch > meta_data.tol && meta_data.iter_count < meta_data.max_iter_num) {
//			obj = MMA_Objective(n_var, var, grad);					//calculate objective and gradient
//			MMA_Constraint(n_cons, cons, cons_grad, n_var, var);	//calculate constraints and gradient
//			Write_Substep(meta_data.iter_count);								//output substep from the initial state
//			//// Set outer move limits
//			if (use_elementwise_bound) {
//				for (int i = 0; i < n_var; i++) {
//					var_U[i] = std::min(var_UB[i], var[i] + movlims[i]);
//					var_L[i] = std::max(var_LB[i], var[i] - movlims[i]);
//				}
//			}
//			else {
//				for (int i = 0; i < n_var; i++) {
//					var_U[i] = std::min(var_ub, var[i] + movlim);
//					var_L[i] = std::max(var_lb, var[i] - movlim);
//				}
//			}
//
//			//// Update MMA next step
//			mma->Update(var, grad, cons, cons_grad, var_L, var_U);
//
//			//// Compute inf norm on design change
//			ch = 0.0; for (int i = 0; i < n_var; i++) { ch = std::max(ch, std::abs(var[i] - intmed_var[i])); }
//
//			meta_data.iter_count++;
//		}
//
//		// Deallocate
//		Delete(mma);
//	}
//
//	void Set_Var_LB(const real* _var_lb)
//	{
//		for (int i = 0; i < n_var; i++) {
//			var_LB[i] = _var_lb[i];
//		}
//		use_var_L = true;
//	}
//
//	void Set_Var_UB(const real* _var_ub)
//	{
//		for (int i = 0; i < n_var; i++) {
//			var_UB[i] = _var_ub[i];
//		}
//		use_var_U = true;
//	}
//
//	void Set_Var_LB_And_UB(const real* _var_lb, const real* _var_ub)
//	{
//		Set_Var_LB(_var_lb);
//		Set_Var_UB(_var_ub);
//	}
//
//public:
//	real MMA_Objective(unsigned n_var, const real* var, real* grad) //Objective and gradient
//	{
//		real obj_value = Compute_Objective(var);
//		std::memcpy(&intmed_var[0], var, n_var * sizeof(real)); //must be computed to check convergence
//		Compute_Gradient(var, grad);
//
//		return obj_value;
//	}
//
//	void MMA_Constraint(unsigned n_cons, real* cons, real* cons_grad, unsigned n_var, const real* var) //Constraint and gradient
//	{
//		if (n_cons > 0) {
//			Compute_Constraint(var, cons);
//			Compute_Constraint_Grad(var, cons_grad);
//		}
//	}
//};
