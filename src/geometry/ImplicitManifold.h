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

	template<int d>
	class ImplicitManifold {
		Typedef_VectorD(d);
	public:
		virtual real Phi(const VectorD pos)const = 0;
		virtual bool Inside(const VectorD pos)const { return Phi(pos) < (real)0; }
		virtual VectorD Normal(const VectorD pos) const = 0;//MUST BE NORMALIZED
	};

	template<int d>
	class ImplicitManifoldPtr: public ImplicitManifold<d> {
		Typedef_VectorD(d);
	public:
		std::shared_ptr<ImplicitManifold<d>> data_ptr;
		real Phi(const VectorD pos)const { return data_ptr->Phi(pos); }
		VectorD Normal(const VectorD pos) const { return data_ptr->Normal(pos); }
	};

	template<int d, class Shape>
	class ImplicitManifoldShape: public ImplicitManifoldPtr<d> {
	public:
		template<class... Args>
		ImplicitManifoldShape(const Args...args) {
			this->data_ptr = std::make_shared<Shape>(args...);
		}
	};

	//The sub-manifolds can't intersect
	template<int d>
	class ImplicitUnion : public ImplicitManifold<d> {
		Typedef_VectorD(d);
	public:
		Array<ImplicitManifoldPtr<d>> manifolds;
		ImplicitUnion(const ImplicitManifoldPtr<d> manifold1, const ImplicitManifoldPtr<d> manifold2) {
			manifolds.push_back(manifold1);
			manifolds.push_back(manifold2);
		}
		std::tuple<real, ImplicitManifoldPtr<d>> Closest_Manifold(const VectorD pos)const {
			//return the one with the minimum abs phi
			real min_phi = std::numeric_limits<real>::max();
			ImplicitManifoldPtr<d> min_manifold;
			for (auto manifold : manifolds) {
				real phi = manifold.Phi(pos);
				if (std::fabs(phi) < std::fabs(min_phi)) {
					min_phi = phi;
					min_manifold = manifold;
				}
			}
			return std::make_tuple(min_phi, min_manifold);
		}
		real Phi(const VectorD pos)const {
			auto [phi, manifold] = Closest_Manifold(pos);
			return phi;
		}
		VectorD Normal(const VectorD pos)const {
			auto [phi, manifold] = Closest_Manifold(pos);
			return manifold.Normal(pos);
		}
	};

	template<int d> 
	class SphereShape : public ImplicitManifold<d>
	{
		Typedef_VectorD(d);
	public:
		VectorD center = VectorD::Zero();
		real radius = (real)1;
		SphereShape(const VectorD _center, const real _radius) :center(_center), radius(_radius) {}
		//SphereShape<d>& operator=(const Sphere<d>& copy) { center = copy.center; radius = copy.radius; return *this; }
		//SphereShape(const Sphere<d>& copy) { *this = copy; }

		real Phi(const VectorD pos) const { return (pos - center).norm() - radius; }
		VectorD Normal(const VectorD pos) const { return (pos - center).normalized(); }
	};

	template<int d>
	class PlaneShape : public ImplicitManifold<d>
	{
		Typedef_VectorD(d);
	public:
		VectorD p;
		VectorD n;
		real b;

		PlaneShape(const VectorD _p, const VectorD _n = VectorD::Unit(1)) :p(_p), n(_n) { n.normalize(); b = n.dot(p); }
		//SphereShape<d>& operator=(const SphereShape<d>& copy) { n = copy.n; p = copy.p; b = copy.b; return *this; }
		//SphereShape(const SphereShape<d>& copy) { *this = copy; }

		bool Inside(const VectorD pos) const { return n.dot(pos) - b < (real)0; }
		real Phi(const VectorD pos) const { return (n.dot(pos) - b); }
		VectorD Normal(const VectorD pos) const { return n; }
	};

    template<int d> class BoxShape : public ImplicitManifold<d>
    {
        Typedef_VectorD(d);
    public:
        VectorD min_corner, max_corner;

		BoxShape(const VectorD _min = VectorD::Zero(), const VectorD _max = VectorD::Zero()) :min_corner(_min), max_corner(_max) {}
		//BoxShape(const VectorD& center, const real side_length) {
  //          VectorD offset = VectorD::Ones() * side_length * 0.5;
  //          min_corner = center - offset;
  //          max_corner = center + offset;
  //      }
		//BoxShape<d>& operator=(const Box<d>& copy) { min_corner = copy.min_corner; max_corner = copy.max_corner; return *this; }
		//BoxShape(const Box<d>& copy) { *this = copy; }

        bool Inside(const VectorD pos) const { return ArrayFunc::All_Greater_Equal(pos, min_corner) && ArrayFunc::All_Less_Equal(pos, max_corner); }
        real Phi(const VectorD pos) const
        {
            VectorD phi = (pos - (min_corner + max_corner) * 0.5).cwiseAbs() - (real).5 * Edge_Lengths(); VectorD zero = VectorD::Zero();
            if (!ArrayFunc::All_Less_Equal(phi, zero)) { return (phi.cwiseMax(zero)).norm(); }return phi.maxCoeff();
        }
        VectorD Normal(const VectorD pos)const { return Wall_Normal(pos); }
        VectorD Edge_Lengths() const { return max_corner - min_corner; }
        //VectorD Center() const { return (real).5 * (min_corner + max_corner); }
        //Box<d> Enlarged(const Box<d>& box2) const { return Box<d>(MathFunc::Cwise_Min(min_corner, box2.min_corner), MathFunc::Cwise_Max(max_corner, box2.max_corner)); }
        //Box<d> Enlarged(const VectorD& length) const { return Box<d>(min_corner - length, max_corner + length); }
        //Box<d> Rescaled(const real factor) const { VectorD length = Edge_Lengths(); return Box<d>(min_corner - length * factor * (real).5, max_corner + length * factor * (real).5); }
        VectorD Wall_Normal(const VectorD pos) const
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
	class EllipsoidShape : public ImplicitManifold<d>
	{
		Typedef_VectorD(d);
	public:
		VectorD center = VectorD::Zero();
		VectorD radius = MathFunc::V<d>(1., 0.5, 1.);
		EllipsoidShape(const VectorD _center, const VectorD _radius) :center(_center), radius(_radius) {}

		real Phi(const VectorD pos) const { const VectorD diff = pos - center; real k1 = MathFunc::Cwise_Divide(diff, radius).norm(); real k2 = (diff / radius.squaredNorm()).norm(); return k1 * (k1 - 1.0) / k2; }
		VectorD Normal(const VectorD pos) const { return (pos - center).normalized(); }	// TODO: fix ellipsoid normal
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

	template<int d> using Sphere = ImplicitManifoldShape<d, SphereShape<d>>;
	template<int d> using Plane = ImplicitManifoldShape<d, PlaneShape<d>>;
	template<int d> using Box = ImplicitManifoldShape<d, BoxShape<d>>;
	template<int d> using Ellipsoid = ImplicitManifoldShape<d, EllipsoidShape<d>>;
}
