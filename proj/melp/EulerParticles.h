//////////////////////////////////////////////////////////////////////////
// Basic grid data representation
// Copyright (c) (2022-), Bo Zhu, Mengdi Wang, Yitong Deng
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "NAParticles.h"
#include "ParticleFrames.h"
#include "Neighbors.h"

namespace Meso {
	template<int d>
	class EulerParticles : public NAParticles<d> {
		Typedef_VectorD(d); Typedef_MatrixD(d);
	public:
		EulerParticles() : NAParticles<d>()  {}

		Setup_Attribute(rho, real, 1.0);
		Setup_Attribute(eta, real, 1.0);
		Setup_Attribute(nden, real, 1.0);
		Setup_Attribute(a, real, 1.0);
		Setup_Attribute(H, real, 0.0); //mean curvature
		Setup_Attribute(u, VectorD, VectorD::Zero());
		Setup_Attribute(E, MatrixD, MatrixD::Identity());
		Setup_Attribute(g, MatrixT, MatrixT::Identity());

		Array<Neighbors> nbs_info;

		real dx = 1.;

		// normal vector i
		VectorD Normal(int i) {
			return E(i).col(d-1);
		}

		// metric tensor of i
		MatrixT Metric(int i) {
			return g(i).inverse();
		}

		// search radius for i
		real Radius(int i) {
			return 4 * dx;
		}

		void Get_Nbs_Pos(Array<int>& nbs, Array<VectorD>& nbs_pos) {
			nbs_pos.resize(nbs.size());
			for (int k = 0; k < nbs.size(); k++) {
				nbs_pos[k] = x(nbs[k]);
			}
		}

		void Update_Local_Frames(void) {
#pragma omp parallel for
			for (int i = 0; i < Size(); i++) {
				Array<int> nbs; Array<VectorD> nbs_pos; 
				nbs_searcher->Find_Nbs(x(i), Radius(i), nbs); 
				Get_Nbs_Pos(nbs, nbs_pos);
				Set_Local_Frame_PCA<d>(Radius(i), x(i), nbs_pos, E(i));
			}
		}

		void Orient_Normals(void) {
			Orient_Normals_COM<d>(xRef(), ERef());
		}
	};
}