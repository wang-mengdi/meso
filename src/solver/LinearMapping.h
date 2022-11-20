#pragma once

#include "Common.h"
#include "AuxFunc.h"

namespace Meso {
	template<class T>
	class LinearMapping
	{
	public:
		virtual int XDoF() const = 0;//former xdof, size of input

		virtual int YDoF() const = 0;//number of rows, size of output

		//input p, get Ap
		//Ap must be allocated, but may be set to arbitrary values
		virtual void Apply(ArrayDv<T>& Ap, const ArrayDv<T>& p) = 0;

		bool Empty(void)const { return XDoF() * YDoF() == 0; }

		void Residual(ArrayDv<T>& res, const ArrayDv<T>& x, const ArrayDv<T> &b) {
			res.resize(b.size());
			//b-Ax
			Apply(res, x);
			//ArrayFunc::Binary_Transform(res, b, [=]__device__(T a, T b) { return b - a; }, res);
			T* res_ptr = thrust::raw_pointer_cast(res.data());
			const T* b_ptr = thrust::raw_pointer_cast(b.data());
			auto f = [=]__device__(T & a, const T & b) { a = b - a; };
			GPUFunc::Cwise_Mapping_Wrapper(res_ptr, b_ptr, f, res.size());
		}

		//check if Ap and p has enough space
		template<typename ...Args>
		void Memory_Check(const ArrayDv<T>& Ap, const ArrayDv<T>& p, const Args&...args) const {
			Assert(p.size() == XDoF() && Ap.size() == YDoF(), args...);
		}
	};

}