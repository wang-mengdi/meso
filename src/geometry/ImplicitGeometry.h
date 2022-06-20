//////////////////////////////////////////////////////////////////////////
// Geometry primitives
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include <limits>
#include <iostream>
#include "Common.h"
#include "AuxFunc.h"
#include "Constants.h"

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

    template<int d> class Box : public ImplicitGeometry<d>
    {
        Typedef_VectorD(d);
    public:
        VectorD min_corner, max_corner;

        Box(const VectorD& _min = VectorD::Zero(), const VectorD& _max = VectorD::Zero()) :min_corner(_min), max_corner(_max) {}
        Box(const VectorD& center, const real side_length) {
            VectorD offset = VectorD::Ones() * side_length * 0.5;
            min_corner = center - offset;
            max_corner = center + offset;
        }
        Box<d>& operator=(const Box<d>& copy) { min_corner = copy.min_corner; max_corner = copy.max_corner; return *this; }
        Box(const Box<d>& copy) { *this = copy; }

        virtual bool Inside(const VectorD& pos) const { return ArrayFunc::All_Greater_Equal(pos, min_corner) && ArrayFunc::All_Less_Equal(pos, max_corner); }
        virtual real Phi(const VectorD& pos) const
        {
            VectorD phi = (pos - Center()).cwiseAbs() - (real).5 * Edge_Lengths(); VectorD zero = VectorD::Zero();
            if (!ArrayFunc::All_Less_Equal(phi, zero)) { return (phi.cwiseMax(zero)).norm(); }return phi.maxCoeff();
        }
        virtual VectorD Normal(const VectorD& pos)const { return Wall_Normal(pos); }
        VectorD Edge_Lengths() const { return max_corner - min_corner; }
        //VectorD Center() const { return (real).5 * (min_corner + max_corner); }
        //Box<d> Enlarged(const Box<d>& box2) const { return Box<d>(MathFunc::Cwise_Min(min_corner, box2.min_corner), MathFunc::Cwise_Max(max_corner, box2.max_corner)); }
        //Box<d> Enlarged(const VectorD& length) const { return Box<d>(min_corner - length, max_corner + length); }
        //Box<d> Rescaled(const real factor) const { VectorD length = Edge_Lengths(); return Box<d>(min_corner - length * factor * (real).5, max_corner + length * factor * (real).5); }
        VectorD Wall_Normal(const VectorD& pos) const
        {
            VectorD normal = VectorD::Zero();
            for (int i = 0; i < d; i++) {
                if (pos[i] < min_corner[i])normal[i] = min_corner[i] - pos[i];
                else if (pos[i] > max_corner[i])normal[i] = max_corner[i] - pos[i];
            }
            if (normal != VectorD::Zero()) { return normal.normalized(); }
            return VectorD::Zero();
        }
        //static Box<d> Infi_Min() { const real fmax = std::numeric_limits<real>::max(); return Box<d>(Vector<real, d>::Ones() * fmax, Vector<real, d>::Ones() * (real)-fmax); }
    };

	template<int d> 
	class Plane : public ImplicitGeometry<d>
	{
		Typedef_VectorD(d);
	public:
		VectorD n;
		VectorD p;
		real b;

		Plane(const VectorD _n, const VectorD _p) :n(_n), p(_p) { n.normalize(); b = n.dot(p); }
		Plane<d>& operator=(const Plane<d>& copy) { n = copy.n; p = copy.p; b = copy.b; return *this; }
		Plane(const Plane<d>& copy) { *this = copy; }

		virtual bool Inside(const VectorD& pos) const { return n.dot(pos) - b < (real)0; }
		virtual real Phi(const VectorD& pos) const { return (n.dot(pos) - b); }
		virtual VectorD Normal(const VectorD& pos) const { return n; }
	};

	//template<int d>
	//class GeometryUnion : public ImplicitGeometry<d> {
	//public:
	//	GeometryUnion(const ImplicitGeometry<d>& a, const ImplicitGeometry<d>& b) {

	//	}
	//	Array<std::shared_ptr<ImplicitGeomtry<d>>> geoms;
	//	virtual real Phi(const VectorD& pos) const {
	//		return (pos - center).norm() - radius;
	//	}
	//	//virtual VectorD Normal(const VectorD& pos) const { return (pos - center).normalized(); }
	//};
}