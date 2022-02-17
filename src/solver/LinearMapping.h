#pragma once

#include "Common.h"

namespace Meso {

	template<class T>
	class LinearMapping
	{
	public:
		virtual int cols() const = 0;//former xdof

		virtual int Y_DoF() const = 0;//number of rows

		//input p, get Ap
		virtual void Apply(ArrayDv<T>& Ap, const ArrayDv<T>& p) = 0;
	};

}