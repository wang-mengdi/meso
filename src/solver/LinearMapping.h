#pragma once

#include "Common.h"

template<class T>
class LinearMapping
{
public:
	virtual int xDoF() const = 0;

	virtual int yDoF() const = 0;

	//input p, get Ap
	virtual void applyMapping(ArrayD& Ap, const ArrayD& p) = 0;
};