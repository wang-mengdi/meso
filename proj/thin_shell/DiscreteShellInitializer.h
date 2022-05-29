//#####################################################################
// Soft Body Nonlinear Thin Shell Driver
// Copyright (c) (2021-), Fan Feng, fan.feng.gr@dartmouth.edu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//#####################################################################
#pragma once
#include <fstream>
#include "AuxFunc.h"
#include "MeshFunc.h"
#include "DiscreteShell.h"

template<int d> class DiscreteShellInitializer
{Typedef_VectorD(d);Typedef_MatrixD(d);
public:
	SurfaceMesh<d> mesh;

	virtual void Apply(json& j, DiscreteShell<d>& thin_shell_simulator)
	{
		int test = Json::Value(j, "test", 0);

		std::string mode = Json::Value(j, "mode", std::string("implicit"));
		if (mode == std::string("explicit")) {thin_shell_simulator.advance_mode = AdvanceMode::Explicit;}
		else if (mode == std::string("implicit")) {thin_shell_simulator.advance_mode = AdvanceMode::Implicit;}
		else if(mode== std::string("quasistatic")){thin_shell_simulator.advance_mode = AdvanceMode::Quasistatic;}
		else{Error("No such advance mode available");}

		switch (test) {
		case 0:Case_0(j, thin_shell_simulator); break;
		case 1:Case_1(j, thin_shell_simulator); break;
		default:Assert(false, "test {} not exist", test); break;
		}
	}
	
	void Case_0(json& j, DiscreteShell<d>& thin_shell_simulator) {
		int scale = Json::Value(j, "scale", 32);
		real length = (real)0.1; int w = scale; real step = length / (real)w;
		MeshFunc::Initialize_Herring_Bone_Mesh(w, w, step, &mesh, 0, 2);
		thin_shell_simulator.Initialize(mesh);

		for (int i = 0; i < w; i++) {
			thin_shell_simulator.Set_Fixed(i);
		}
		thin_shell_simulator.use_body_force = true;

		thin_shell_simulator.damping = Json::Value(j, "damping", (real)0.01);
		thin_shell_simulator.thickness = Json::Value(j, "thickness", (real)0.01);
		real youngs = Json::Value(j, "youngs", (real)100);
		real poisson = Json::Value(j, "poisson", (real)0.35);
		real density = Json::Value(j, "density", (real)1e3);

		thin_shell_simulator.materials.clear();
		thin_shell_simulator.Add_Material((real)youngs, (real)poisson);
		thin_shell_simulator.density = density;
		ArrayFunc::Fill(thin_shell_simulator.material_id, 0);
		thin_shell_simulator.Initialize_Material();
	}

	void Case_1(json& j, DiscreteShell<d>& thin_shell_simulator) {
		int scale = Json::Value(j, "scale", 32);
		real length = (real)0.1; int w = scale; real step = length / (real)w;
		MeshFunc::Initialize_Herring_Bone_Mesh(w, w, step, &mesh, 0, 2);
		thin_shell_simulator.Initialize(mesh);

		real strength = Json::Value(j, "strength", (real)0.1);
		for (int i = 0; i < w; i++) {
			thin_shell_simulator.Set_Fixed(i);
		}
		for (int i = thin_shell_simulator.particles.Size() - w; i < thin_shell_simulator.particles.Size(); i++) { thin_shell_simulator.Set_Force(i, -strength*VectorD::Unit(1)); }

		thin_shell_simulator.damping = Json::Value(j, "damping", (real)0.01);
		thin_shell_simulator.thickness = Json::Value(j, "thickness", (real)0.01);
		real youngs = Json::Value(j, "youngs", (real)100);
		real poisson = Json::Value(j, "poisson", (real)0.35);
		real density = Json::Value(j, "density", (real)1e3);

		thin_shell_simulator.materials.clear();
		thin_shell_simulator.Add_Material((real)youngs, (real)poisson);
		thin_shell_simulator.density = density;
		ArrayFunc::Fill(thin_shell_simulator.material_id, 0);
		thin_shell_simulator.Initialize_Material();
	}
};