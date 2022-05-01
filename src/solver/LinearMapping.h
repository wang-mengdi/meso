#pragma once

#include "Common.h"
#include "AuxFunc.h"

namespace Meso {
	template<class T>
	class LinearMapping
	{
	public:
		virtual int XDoF() const = 0;//former xdof

		virtual int YDoF() const = 0;//number of rows

		//input p, get Ap
		//Ap must be allocated, but may be set to arbitrary values
		virtual void Apply(ArrayDv<T>& Ap, const ArrayDv<T>& p) = 0;

		void Residual(ArrayDv<T>& res, const ArrayDv<T>& x, const ArrayDv<T> &b) {
			//b-Ax
			Apply(res, x);
			ArrayFunc::Binary_Transform(res, b, [=]__device__(T a, T b) { return b - a; }, res);
		}

		//check if Ap and p has enough space
		template<typename ...Args>
		void Memory_Check(const ArrayDv<T>& Ap, const ArrayDv<T>& p, const Args&...args) const {
			Assert(p.size() == XDoF() && Ap.size() == YDoF(), args...);
		}
	};

}