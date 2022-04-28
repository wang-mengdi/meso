//////////////////////////////////////////////////////////////////////////
// Geometry primitives
// Copyright (c) (2018-), Bo Zhu, Wanxin Hu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __GeometryPrimitives_h__
#define __GeometryPrimitives_h__
#include <limits>
#include <iostream>
#include "Common.h"
#include "AuxFunc.h"
#include "Constants.h"

namespace Meso {

    class ImplicitGeometry {
        Typedef_VectorDii(d);
    public:
        ImplicitGeometry() {};
        virtual real Phi(const VectorD& pos) const = 0;
        virtual bool Inside(const VectorD& pos) const {
            return Phi(pos) < (real)0;
        }
        virtual VectorD Normal(const VectorD& pos) const = 0;
    }

    template<int d> class Sphere : public ImplicitGeometry<d> {
        Typedef_VectorDii(d);
    public:
        VectorD center = VectorD::Zero();
        real radius = (real)1;

        Sphere(const VectorD& _radius, const VectorD& _center) :center(_center), radius(_radius) {}
        Sphere<d>& operator=(const Sphere<d>& copy) { center = copy.centerl; radius = copy.radius; return *this; }
        Sphere(const Sphere<d> copy) { *this = copy; }

        virtual real Phi(const VectorD& pos) const {
            return (pos - center).norm() - radius;
        }
        virtual VectorD Normal(const VectorD& pos) const {
            return (pos - center).normalized();
        }
    };
}
#endif