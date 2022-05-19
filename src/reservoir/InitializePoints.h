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
	using namespace CommonConstants;

	void Initialize_Lattice_Points(const Vector3& origin, const int nx, const int ny, const Vector3 dx, const Vector3 dy, Points& pts, Array<Vector3>& pos) {
		pts.Resize((nx + 1) * (ny + 1));
		int idx = 0;
#pragma omp parallel for
		for (int i = 0; i <= nx; i++) {
			for (int j = 0; j <= ny; j++) {
				Vector3 curr_pos = origin + i * dx + j * dy;
				pos[idx] = curr_pos;
				idx++;
			}
		}
	}


	// algorithm: https://www.cmu.edu/biolphys/deserno/pdf/sphere_equi.pdf
	// note that the num_pts you put in is only an average, it is not necessarily the num of pts you get
	real Initialize_Sphere_Points_Regular(const Vector3& origin, const real R, const int num_pts, Points& pts, Array<Vector3>& pos) {
		Array<real> dxs;
		Array<Vector3> tmp_pos;
		real a = 4 * pi / num_pts;
		real d = sqrt(a);
		int M_theta = (int)round(pi/d);
		real d_theta = pi / M_theta;
		real d_phi = a / d_theta;
#pragma omp parallel for
		for (int m = 0; m < M_theta; m++) {
			real theta = pi * ((real)m + 0.5) / M_theta;
			int M_phi = (int)round(2 * pi * sin(theta) / d_phi);
			for (int n = 0; n < M_phi; n++) {
				real phi = 2 * pi * n / M_phi;
				real x = sin(theta) * cos(phi);
				real y = sin(theta) * sin(phi);
				real z = cos(theta);
				Vector3 p(R*x, R*y, R*z);
				p += origin;
				tmp_pos.push_back(p);
			}
			dxs.push_back(0.5 * (d_phi + d_theta));
		}
		real avg_dx = ArrayFunc::Mean<real>(dxs);

		int actual_num_pts = tmp_pos.size();
		pts.Resize(actual_num_pts);
#pragma omp parallel for
		for (int i = 0; i <= actual_num_pts; i++) {
			pos[i] = tmp_pos[i];
		}
		return R * avg_dx;
	}
}