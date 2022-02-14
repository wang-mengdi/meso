#pragma once

#include "Common.h"

template<class T>
class LinearMapping
{
public:
	virtual int xDoF() const = 0;//number of cols

	virtual int yDoF() const = 0;//number of rows

	//input p, get Ap
	virtual void applyMapping(ArrayDv<T>& Ap, const ArrayDv<T>& p) = 0;
};