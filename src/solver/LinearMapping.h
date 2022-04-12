#pragma once

#include "Common.h"
#include "AuxFunc.h"

namespace Meso {
	template<class T>
	class LinearMapping
	{
	public:
		virtual int XDof() const = 0;//former xdof

		virtual int YDof() const = 0;//number of rows

		//input p, get Ap
		virtual void Apply(ArrayDv<T>& Ap, const ArrayDv<T>& p) = 0;

		virtual void Residual(ArrayDv<T>& res, const ArrayDv<T>& x, const ArrayDv<T> &b) {
			//b-Ax
			Apply(res, x);
			ArrayFunc::Binary_Transform(res, b, [=]__device__(T a, T b) { return b - a; }, res);
		}

		virtual bool Size_Match(const ArrayDv<T>& Ap, const ArrayDv<T>& p) const {
			return p.size() == XDof() && Ap.size() == YDof();
		}
	};

}