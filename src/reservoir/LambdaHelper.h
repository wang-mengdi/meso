//////////////////////////////////////////////////////////////////////////
// Functional programming helpers
// Copyright (c) (2022-), Mengdi
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "Common.h"

namespace Meso {

	//I for index
	//C for cell coord (VectorDi)
	//F for face coord (VectorDi)

	template<class T> using IFunc = std::function<T(const int i)>;
	template<class T, int d> using CFunc = std::function<T(const Vector<int, d> cell)>;
	template<class T, int d> using IFFunc = std::function<T>(const int axis, const Vector<int, d> face);

	template<class Fint>
    class Functor_Wrapper
    {

        const float a;

        saxpy_functor(float _a) : a(_a) {
        }

        __host__ __device__
            float operator()(const float& x, const float& y) const {

            return a * x + y;
        }
    };

}