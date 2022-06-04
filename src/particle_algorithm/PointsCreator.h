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
	// initialize planar points in 3D
	void Initialize_Lattice_Points(const Vector3& origin, const int nx, const int ny, const Vector3 dx, const Vector3 dy, Points& pts, Array<Vector3>& pos);
	// 
	real Initialize_Box_Points_2D(const Vector2& origin, const Vector2i& size, const Vector2& _dx, Points& pts, Array<Vector2>& pos, bool keep_existing = false);
	real Initialize_Box_Rim_Points_2D(const Vector2& origin, const Vector2i& pad_size, const Vector2i& int_size, const Vector2& _dx, Points& pts, Array<Vector2>& pos, bool keep_existing = false);
	real Initialize_Box_Points_3D(const Vector3& origin, const Vector3i& size, const Vector3& _dx, Points& pts, Array<Vector3>& pos, bool keep_existing = false);
	real Initialize_Box_Rim_Points_3D(const Vector3& origin, const Vector3i& pad_size, const Vector3i& int_size, const Vector3& _dx, Points& pts, Array<Vector3>& pos, bool keep_existing = false);
	//void Initialize_Box_Points_3D()
	real Initialize_Sphere_Points_Regular(const Vector3& origin, const real R, const int num_pts, Points& pts, Array<Vector3>& pos);
}