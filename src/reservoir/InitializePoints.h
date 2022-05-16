//////////////////////////////////////////////////////////////////////////
// Point initialization Functions
// Copyright (c) (2018-), Bo Zhu, Mengdi Wang, Yitong Deng
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "Common.h"
#include "AuxFunc.h"
#include "Points.h"

namespace Meso {
	void Initialize_Lattice_Points(const Vector3& origin, const int nx, const int ny, const real dx, const real dy, Points& pts, std::string name) {
		try {
			Array<Vector3>& pos = pts.Get_Attribute<Vector3>("x");
			pts.Resize((nx + 1) * (ny + 1));
			int idx = 0;
			for (int i = 0; i <= nx; i++) {
				for (int j = 0; j <= ny; j++) {
					Vector3 curr_pos = origin + Vector3::Unit(0) * i * dx + Vector3::Unit(1) * j * dy;
					pos[idx] = curr_pos;
					idx++;
				}
			}
		}
		catch (const std::exception& e) {
			Error("Error: Failed to extract positions with: {} in InitializePoints.h", name);
		}
	}
}