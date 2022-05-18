//////////////////////////////////////////////////////////////////////////
// Test mesh structure
// Copyright (c) (2022-), Bo Zhu, Yunquan Gu
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once

#include "Mesh.h"
#include "IOFunc.h"

namespace Meso {
	template<class T>
	void Test_Mesh_Loader_Multiple() {
		Array < std::shared_ptr<T>> meshes;
		OBJFunc::Read_Meshes<T>("../../../../../proj/_tests_reservoir/assets/CornellBox-Sphere.obj", meshes);
		OBJFunc::Write_Meshes<T>("./copy-mesh.obj", meshes);
		
		Array < std::shared_ptr<T>> copy_meshes;
		OBJFunc::Read_Meshes<T>("./copy-mesh.obj", copy_meshes);
		Assert(meshes.size() == copy_meshes.size(), "Test_Mesh_Loader Failed: Reload mesh size doesn't match.");
		for (size_t i = 0; i < meshes.size(); i++){
			auto& mesh = meshes[i];
			auto& c_mesh = copy_meshes[i];
			Assert(mesh->Faces().size() == c_mesh->Faces().size(), "Test_Mesh_Loader Failed: faces size doesn't match.", i);
			Assert(mesh->Vertices().size() == c_mesh->Vertices().size(), "Test_Mesh_Loader Failed: vertice size doesn't match.", i);
			for (size_t j = 0; j < mesh->Faces().size(); j++) Assert(mesh->faces[j] == c_mesh->faces[j], "Test_Mesh_Loader Failed: faces value doesn't match.");
			for (size_t j = 0; j < mesh->Vertices().size(); j++) Assert(mesh->Vertices()[j] == c_mesh->Vertices()[j], "Test_Mesh_Loader Failed: vertex value doesn't match.");
		}
		Pass("Test_Mesh_Loader[Multiple] Passed!");
	}

	template<class T>
	void Test_Mesh_Loader_Single() {
		auto mesh = std::make_shared<T>();
		OBJFunc::Read_Mesh("../../../../../proj/_tests_reservoir/assets/CornellBox-Single.obj", mesh);
		OBJFunc::Write_Mesh<T>("./copy_mesh_single.obj", mesh);

		auto copy_mesh = std::make_shared<T>();
		OBJFunc::Read_Mesh("./copy_mesh_single.obj", copy_mesh);
		Assert(mesh->Faces().size() == copy_mesh->Faces().size(), "Test_Mesh_Loader Failed: faces size doesn't match.");
		Assert(mesh->Vertices().size() == copy_mesh->Vertices().size(), "Test_Mesh_Loader Failed: vertice size doesn't match.");
		for (size_t j = 0; j < mesh->Faces().size(); j++) Assert(mesh->faces[j] == copy_mesh->faces[j], "Test_Mesh_Loader Failed: faces value doesn't match.");
		for (size_t j = 0; j < mesh->Vertices().size(); j++) Assert(mesh->Vertices()[j] == copy_mesh->Vertices()[j], "Test_Mesh_Loader Failed: vertex value doesn't match.");
		Pass("Test_Mesh_Loader[Single] Passed!");
	}
}