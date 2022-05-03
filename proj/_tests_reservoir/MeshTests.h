//////////////////////////////////////////////////////////////////////////
// Test mesh structure
// Copyright (c) (2022-), Bo Zhu, Yunquan Gu
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once

#include "Mesh.h"
#include "IOFunc.h"

namespace Meso {
	template<int d> 
	void Test_Mesh_Loader() {
		auto meshes = OBJFunc::Read_Mesh_From_Obj_File("../../../../../proj/_tests_reservoir/assets/CornellBox-Sphere.obj");
		//Info("Test_Mesh_Loader loaded CornellBox-Sphere.obj");
		OBJFunc::Write_Mesh_To_Obj_File("./copy-mesh.obj", meshes);
		//Info("Test_Mesh_Loader writes copy-mesh.obj");
		auto copy_meshes = OBJFunc::Read_Mesh_From_Obj_File("./copy-mesh.obj");
		//Info("Test_Mesh_Loader loaded copy-mesh.obj");
		Assert(meshes.size() == copy_meshes.size(), "Test_Mesh_Loader Failed: Reload mesh size doesn't match.");
		for (size_t i = 0; i < meshes.size(); i++){
			auto& mesh = meshes[i];
			auto& c_mesh = copy_meshes[i];
			Assert(mesh->Faces().size() == c_mesh->Faces().size(), "Test_Mesh_Loader Failed: faces size doesn't match.", i);
			Assert(mesh->Vertices().size() == c_mesh->Vertices().size(), "Test_Mesh_Loader Failed: vertice size doesn't match.", i);
			//Info("Test_Mesh_Loader: check {}-th mesh with {} faces and {} vertices.", i, mesh->Faces().size(), mesh->Vertices().size());
			for (size_t j = 0; j < mesh->Faces().size(); j++) Assert(mesh->faces[j] == c_mesh->faces[j], "Test_Mesh_Loader Failed: faces value doesn't match.");
			for (size_t j = 0; j < mesh->Vertices().size(); j++) Assert(mesh->Vertices()[j] == c_mesh->Vertices()[j], "Test_Mesh_Loader Failed: vertex value doesn't match.");
		}
		Pass("Test_Mesh_Loader Passed!");
	}
}