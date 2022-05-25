//////////////////////////////////////////////////////////////////////////
// Geometry primitives
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include <limits>
#include <iostream>
#include "Common.h"

namespace Meso {
	////////////////////////////////////////////////////////////////////////
	//Base class
	////////////////////////////////////////////////////////////////////////
	template<int d> class ImplicitGeometry
	{
		Typedef_VectorD(d);
	public:
		ImplicitGeometry() {}
		virtual real Phi(const VectorD& pos) const = 0;
		virtual bool Inside(const VectorD& pos) const { return Phi(pos) < (real)0; }
		virtual VectorD Normal(const VectorD& pos) const = 0;//MUST BE NORMALIZED
	};

	template<int d> class Sphere : public ImplicitGeometry<d>
	{
		Typedef_VectorD(d);
	public:
		VectorD center = VectorD::Zero();
		real radius = (real)1;
		Sphere(const VectorD& _center, const real _radius) :center(_center), radius(_radius) {}
		Sphere<d>& operator=(const Sphere<d>& copy) { center = copy.center; radius = copy.radius; return *this; }
		Sphere(const Sphere<d>& copy) { *this = copy; }

		virtual real Phi(const VectorD& pos) const { return (pos - center).norm() - radius; }
		virtual VectorD Normal(const VectorD& pos) const { return (pos - center).normalized(); }
	};
}