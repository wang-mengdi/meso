//#####################################################################
// Soft Body Nonlinear Thin Shell Driver
// Copyright (c) (2021-), Fan Feng, fan.feng.gr@dartmouth.edu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//#####################################################################
#pragma once
#include <fstream>
#include "AuxFunc.h"
#include "MeshFunc.h"
#include "SoftBodyNonlinearFemThinShell.h"

template<int d> class SoftBodyNonlinearThinShellInitializer
{Typedef_VectorD(d);Typedef_MatrixD(d);
public:
	SurfaceMesh<d> mesh;

	virtual void Apply(json& j, SoftBodyNonlinearFemThinShell<d>& thin_shell_simulator)
	{
		int test = Json::Value(j, "test", 0);
		int scale = Json::Value(j, "scale", 32);
		Initialize_Mesh(test,scale, thin_shell_simulator);
		Initialize_Boundary_Conditions(test,scale, thin_shell_simulator);
		Initialize_Materials(test, thin_shell_simulator);
	}

	void Initialize_Mesh(int test,int scale, SoftBodyNonlinearFemThinShell<d>& thin_shell_simulator)
	{
		if constexpr (d == 2) {
			switch (test) {
			case 1: {
				MeshFunc::Initialize_Segment_Mesh<d>(VectorD::Zero(), VectorD::Unit(0), scale, &mesh,false,true);
			}break;
			}
		}

		if constexpr (d == 3) {
			switch (test) {
			case 1: {
				real length = (real)1; int w = scale; real step = length / (real)w;
				MeshFunc::Initialize_Herring_Bone_Mesh(w, w, step, &mesh, 0, 2);
			}break;
			}
		}

		thin_shell_simulator.Initialize(mesh);
	}

	void Initialize_Boundary_Conditions(int test,int scale, SoftBodyNonlinearFemThinShell<d>& thin_shell_simulator)
	{
		switch(test){
		case 1:{	////beam with one end fixed under gravity
			if constexpr (d == 2) {
				thin_shell_simulator.Set_Fixed(0);
			}
			if constexpr (d == 3) {
				real length = (real)1; int w = scale; real step = length / (real)w;

				for (int i = 0; i < w; i++) {
					thin_shell_simulator.Set_Fixed(i);
				}
			}
			thin_shell_simulator.use_body_force=true;
		}break;
		}
	}

	void Initialize_Materials(int test, SoftBodyNonlinearFemThinShell<d>& thin_shell_simulator)
	{
		switch (test) {
			case 1:{
				thin_shell_simulator.materials.clear();
				thin_shell_simulator.Add_Material((real)1e2, (real).35);
				ArrayFunc::Fill(thin_shell_simulator.material_id, 0);
				thin_shell_simulator.Initialize_Material();
			}break;
		}
	}
};