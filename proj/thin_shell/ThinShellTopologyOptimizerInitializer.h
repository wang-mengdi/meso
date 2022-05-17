//////////////////////////////////////////////////////////////////////////
// Topology Optimization for Thin Shell
// Copyright (c) (2021-), Fan Feng
// This file is part of CompleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include <fstream>
#include "Common.h"
#include "Driver.h"
#include "AuxFunc.h"
#include "MeshFunc.h"
#include "ThinShellTopologyOptimizer.h"
#include <sstream>

template<int d> class TopoOptThinShellInitializer
{Typedef_VectorD(d); Typedef_MatrixD(d);

public:
	SurfaceMesh<d> mesh;
	int w;
	int h;
	int test;

	virtual void Apply(json& j, ThinShellTopologyOptimizer<d>& optimizer)
	{
		test = Json::Value(j, "test", 0);
		w = Json::Value(j, "width", 32);
		h = Json::Value(j, "height", 32);
		Initialize_Mesh(optimizer);
		Initialize_Boundary_Conditions(optimizer);
		Initialize_Materials(optimizer,j);
		Initialize_Topo_Opt(j,optimizer);
	}

	void Initialize_Mesh(ThinShellTopologyOptimizer<d>& optimizer)
	{
		SoftBodyNonlinearFemThinShellTopoOpt<d>& thin_shell = optimizer.thin_shell;

		if constexpr (d == 3) {
			switch (test) {
				case 1: 
				case 5:
				case 7:
				case 8:
				case 9:
				{
					real length = (real)1; real step = length / (real)w;
					MeshFunc::Initialize_Herring_Bone_Mesh(w, h, step, &mesh, 0, 2);
				}break;
				case 2: {
					real length = (real)1; real step = length / (real)w;
					MeshFunc::Initialize_Herring_Bone_Mesh(w, h, step, &mesh, 0, 2);
					mesh.Vertices()[0] -= 0.01 * VectorD::Unit(1);
					mesh.Vertices()[w * h] -= 0.01 * VectorD::Unit(1);
					mesh.Vertices()[w-1] += 0.01* VectorD::Unit(1);
					mesh.Vertices()[w*h-w] += 0.01 * VectorD::Unit(1);
				}break;
				case 3: {
					real length = (real)1; real step = length / (real)w;
					MeshFunc::Initialize_Herring_Bone_Mesh(w, h, step, &mesh, 0, 2);
					//for (int i = w; i < w * h - w; i++) { mesh.Vertices()[i] += 0.1 * VectorD::Unit(1); }
				}break;
				case 4: {
					real length = (real)1; real step = length / (real)w;
					MeshFunc::Initialize_Herring_Bone_Mesh(w, h, step, &mesh, 0, 2);
					for (int i = w; i < w * h - w; i++) { mesh.Vertices()[i] += 0.01 * VectorD::Unit(1); }
				}break;
				case 6: {
					real length = (real)1; real step = length / (real)w;
					MeshFunc::Initialize_Herring_Bone_Mesh(w, h, step, &mesh, 0, 2);
					//for (int i = w * h - w; i < w * h; i++) { mesh.Vertices()[i] += 0.1 * VectorD::Unit(1); }
				}break;
			}
		}
		thin_shell.use_exact_hessian = true;
		thin_shell.Initialize(mesh);
	}

	void Initialize_Boundary_Conditions(ThinShellTopologyOptimizer<d>& optimizer)
	{
		SoftBodyNonlinearFemThinShellTopoOpt<d>& thin_shell = optimizer.thin_shell;
		if constexpr (d == 3) {
			switch (test) {
				case 1: {	////beam with one end fixed (two rows) and push on the other side
					real strength = 1e-5;
					for (int i = 0; i < 2 * w; i++) { thin_shell.Set_Fixed(i); }
					//for (int i = thin_shell.particles.Size() - w; i < thin_shell.particles.Size(); i++) { thin_shell.Set_Force(i, -strength *((real)(i-thin_shell.particles.Size() + w)/(real)w)* VectorD::Unit(1) / (real)w); }
					//for (int i = thin_shell.particles.Size() - w; i < thin_shell.particles.Size(); i++) { thin_shell.Set_Force(i, -strength * VectorD::Unit(1)/ (real)w) ; }
					thin_shell.Set_Force(w*h-1, -strength * VectorD::Unit(1));
				}break;
				case 2: {	////pull two ends up and other two ends down, points
					real strength = 0.1;
					thin_shell.Set_Force(0, -strength * VectorD::Unit(1) / (real)w);
					thin_shell.Set_Force(w * w - 1, -strength * VectorD::Unit(1) / (real)w);
					thin_shell.Set_Force(w - 1, strength * VectorD::Unit(1) / (real)w);
					thin_shell.Set_Force(w * w - w, strength * VectorD::Unit(1) / (real)w);
				}break;
				case 3: {	////pull two ends out, edges
					real strength = 0.1;
					for (int i = 0; i < w; i++) { thin_shell.Set_Force(i, -strength * VectorD::Unit(2) / (real)w); }
					for (int i = w*h-w; i < w*h; i++) { thin_shell.Set_Force(i, strength * VectorD::Unit(2) / (real)w); }
					//for (int i = w; i < w * h - w; i++) { thin_shell.Set_Fixed(i); }
				}break;
				case 4: {  ////fix one end and push one end in with dirichlet boundary condition
					real strength = 2;
					for (int i = 0; i < w; i++) { thin_shell.Set_Displacement(i, VectorD::Unit(2)*(real)strength); }
					for (int i = w * h - w; i <w*h; i++) { thin_shell.Set_Displacement(i, -0.1*VectorD::Unit(2) * (real)strength); }
				}break;
				case 5: { //fix one end and pull out the other end
					real strength = 0.1;
					for (int i = 0; i < w; i++) { thin_shell.Set_Fixed(i); }
					for (int i = w * h - w; i < w * h; i++) { thin_shell.Set_Force(i, strength * VectorD::Unit(2) / (real)w); }
				}break;
				case 6: {//fix one end and pull in a tilted direction
					real strength = 0.01;
					for (int i = 0; i < 2 * w; i++) { thin_shell.Set_Fixed(i); }
					for (int i = w * h - w; i < w * h; i++) { thin_shell.Set_Force(i, strength *  VectorD::Unit(1) / (real)w); }
				}break;
				case 7: {//in-plane cantilever beam
					real strength = 0.0001;
					for (int i = 0; i < w; i++) { thin_shell.Set_Fixed(i); }
					thin_shell.Set_Force(w * h - w, strength * -VectorD::Unit(0));
				}break;
				case 8: {//fix four points and add force in the middle
					real strength = 0.1;
					thin_shell.Set_Fixed(0); thin_shell.Set_Fixed(w - 1); thin_shell.Set_Fixed(w * h - w); thin_shell.Set_Fixed(w * h - 1);
					thin_shell.Set_Force(w * (h / 2)+w/2, strength * -VectorD::Unit(1));
				}break;
				case 9: {//Stretch four ends out from center
					real strength = 0.2;
					thin_shell.Set_Fixed(w * (h / 2) + w / 2);
					VectorD middle = mesh.Vertices()[w * (h / 2) + w / 2];
					thin_shell.Set_Force(0, strength * (mesh.Vertices()[0]-middle).normalized());
					thin_shell.Set_Force(w * h-1, strength * (mesh.Vertices()[w * h-1] - middle).normalized());
					thin_shell.Set_Force(w - 1, strength * (mesh.Vertices()[w-1] - middle).normalized());
					thin_shell.Set_Force(w * h - w, strength * (mesh.Vertices()[w * h - w] - middle).normalized());
				}break;
			}
		}
	}

	void Initialize_Materials(ThinShellTopologyOptimizer<d>& optimizer,json& j)
	{
		SoftBodyNonlinearFemThinShellTopoOpt<d>& thin_shell = optimizer.thin_shell;

		switch (test) {
		case 1:
		case 2: 
		case 3: 
		case 4: 
		case 5: 
		case 6: 
		case 7:
		case 8: 
		case 9: {
			thin_shell.materials.clear();
			thin_shell.Add_Material((real)1e2, (real).35); //this can be tuned
			ArrayFunc::Fill(thin_shell.material_id, 0);
			thin_shell.Initialize_Material();
			ArrayFunc::Fill(thin_shell.hs, optimizer.Rho_To_H(optimizer.target_rho));
		}break;
		}
	}

	void Initialize_Topo_Opt(json& j, ThinShellTopologyOptimizer<d>& optimizer)
	{
		optimizer.target_rho = Json::Value(j, "target_rho", 0.3);
		switch (test) {
		case 6: //sphere with bottom on one side
			std::shared_ptr<ImplicitGeometry<d>> target_shape = std::make_shared<Sphere<d>>(Vector3(0.5, (real)0.5, 0.),(real)0.5);
			optimizer.target_shape = target_shape;
		break;
		}
		optimizer.Init();
	}
};