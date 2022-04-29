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
	bool use_quasi_static=false;

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

	//virtual void Advance_One_Time_Step(const real dt, const real time)
	//{
	//	if (use_quasi_static) {
	//		thin_shell.Advance_Quasi_Static();
	//		//Write_Output_Files(1);
	//		exit(0);
	//	}
	//	else {
	//		thin_shell_simulator.Advance(dt, time);
	//	}
	//}

	//virtual void Write_Output_Files(const int frame)
	//{	
	//	Base::Write_Output_Files(frame);

	//	{std::string file_name=frame_dir+(d==2?"/segment_mesh":"/triangle_mesh");
	//	thin_shell.mesh->Write_To_File_3d(file_name);}

	//	{std::string file_name=frame_dir+"/particles";
	//	thin_shell.particles.Write_To_File_3d(file_name);}

	//	{std::string file_name=frame_dir+"/mat";
	//	int n=(int)thin_shell.material_id.size();
	//	Field<real,1> mat;mat.Resize(n);
	//	for(int i=0;i<n;i++){
	//		mat.array[i]=(real)thin_shell.material_id[i];
	//	}
	//	mat.Write_To_File_3d(file_name);}

	//	{std::string file_name = frame_dir + "/psi_D";
	//	Particles<d> psi_D_particles;
	//	for (auto& p : thin_shell.bc.psi_D_values) {
	//		int idx = p.first;
	//		int i = psi_D_particles.Add_Element(); psi_D_particles.X(i) = thin_shell.particles.X(idx);
	//	}
	//	psi_D_particles.Write_To_File_3d(file_name); }

	//	{std::string file_name = frame_dir + "/psi_N";
	//	Particles<d> psi_N_particles;
	//	for (auto& p : thin_shell.bc.forces) {
	//		int idx = p.first;
	//		int i = psi_N_particles.Add_Element(); psi_N_particles.X(i) = thin_shell.particles.X(idx);
	//		psi_N_particles.F(i) = p.second;
	//	}
	//	psi_N_particles.Write_To_File_3d(file_name); }

	//	if (frame == last_frame) {
	//		std::string file_name = output_dir + "/energy.txt";
	//		std::ofstream energy_file(file_name);
	//		for (int count = 0; count < thin_shell.energies_n.size(); count++) {
	//			energy_file << thin_shell.energies_n[count] << "\n";
	//		}
	//	}
	//}
};