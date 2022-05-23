//////////////////////////////////////////////////////////////////////////
// Point initialization Functions
// Copyright (c) (2018-), Bo Zhu, Mengdi Wang, Yitong Deng
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "Common.h"
#include "AuxFunc.h"
#include "Points.h"
#include "IOFunc.h"
#include "Constants.h"

namespace Meso {
	void Initialize_Lattice_Points(const Vector3& origin, const int nx, const int ny, const Vector3 dx, const Vector3 dy, Points& pts, Array<Vector3>& pos);
	real Initialize_Sphere_Points_Regular(const Vector3& origin, const real R, const int num_pts, Points& pts, Array<Vector3>& pos);
}