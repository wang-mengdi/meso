#pragma once

#include "Common.h"

namespace Meso {

	template<class T>
	class LinearMapping
	{
	public:
		virtual int X_DoF() const = 0;//number of cols

		virtual int Y_DoF() const = 0;//number of rows

		//input p, get Ap
		virtual void Apply(ArrayDv<T>& Ap, const ArrayDv<T>& p) = 0;
	};

}