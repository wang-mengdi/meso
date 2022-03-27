//////////////////////////////////////////////////////////////////////////
// Damped Jacobian Smoother
// Copyright (c) (2022-), Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once

#include "LinearMapping.h"
#include "PoissonFunc.h"

namespace Meso {
	template<class T>
	class Smoother {
	public:
		virtual int Dof()const = 0;
		virtual void Apply(ArrayDv<T> &x_new, const ArrayDv<T>& x_old) = 0;
	};

	template<class T>
	class DampedJacobianSmoother : public Smoother<T> {
	public:
		LinearMapping<T>* mapping;
		T omega;
		int dof;
		ArrayDv<T> diag;
		ArrayDv<T> rhs;
		DampedJacobianSmoother() {}
		template<int d> DampedJacobianSmoother(PoissonMapping<T, d>& _mapping, const ArrayDv<T>& _rhs, const real _omega = 2.0 / 3.0) { Init(_mapping, _rhs, _omega); }
		template<int d>
		void Init(PoissonMapping<T, d>& _mapping, const ArrayDv<T> &_rhs, const T _omega = 2.0 / 3.0) {
			mapping = &_mapping;
			omega = _omega;
			dof = mapping->XDof();
			Poisson_Diagonal(diag, _mapping);
			rhs = _rhs;//deep copy
		}
		virtual int Dof()const { return dof; }
		virtual void Apply(ArrayDv<T>& x_new, const ArrayDv<T>& x_old) {
			//Info("x: {}", x_old);
			//Ax
			mapping->Apply(x_new, x_old);
			//Info("Ax: {}", x_new);
			//b-Ax
			ArrayFunc::Binary_Transform(x_new, rhs, [=]__device__(T a, T b) { return b - a; }, x_new);
			//Info("b-Ax: {}", x_new);
			//(b-Ax)/.rhs
			ArrayFunc::Binary_Transform(x_new, diag, [=]__device__(T a, T b) { return a / b; }, x_new);
			//Info("(b-Ax)/rhs: {}", x_new);
			//x+=(b-Ax)/.rhs*.omega
			real _omega = omega;
			ArrayFunc::Binary_Transform(x_new, x_old, [=]__device__(T a, T b) { return b + a * _omega; }, x_new);
		}
	};
}