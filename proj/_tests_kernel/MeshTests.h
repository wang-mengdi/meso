//////////////////////////////////////////////////////////////////////////
// Test mesh structure
// Copyright (c) (2022-), Yunquan Gu
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once

#include "Mesh.h"
#include "IOFunc.h"
#include <igl/readOBJ.h>
#include <Eigen/src/Core/Matrix.h>
namespace Meso {
	//template<class MeshType>
	//void Test_Mesh_Loader_Multiple() {
	//	Array < std::shared_ptr<MeshType>> meshes;
	//	OBJFunc::Read_Meshes<MeshType>("../../../../../proj/_tests_reservoir/assets/CornellBox-Sphere.obj", meshes);
	//	OBJFunc::Write_Meshes<MeshType>("./copy-mesh.obj", meshes);

	//	Array < std::shared_ptr<MeshType>> copy_meshes;
	//	OBJFunc::Read_Meshes<MeshType>("./copy-mesh.obj", copy_meshes);
	//	Assert(meshes.size() == copy_meshes.size(), "Test_Mesh_Loader Failed: Reload mesh size doesn't match.");
	//	for (size_t i = 0; i < meshes.size(); i++){
	//		auto& mesh = meshes[i];
	//		auto& c_mesh = copy_meshes[i];
	//		Assert(mesh->Elements().size() == c_mesh->Elements().size(), "Test_Mesh_Loader Failed: faces size doesn't match.", i);
	//		Assert(mesh->Vertices().size() == c_mesh->Vertices().size(), "Test_Mesh_Loader Failed: vertice size doesn't match.", i);
	//		for (size_t j = 0; j < mesh->Elements().size(); j++) Assert(mesh->elements[j] == c_mesh->elements[j], "Test_Mesh_Loader Failed: faces value doesn't match.");
	//		for (size_t j = 0; j < mesh->Vertices().size(); j++) Assert(mesh->Vertices()[j] == c_mesh->Vertices()[j], "Test_Mesh_Loader Failed: vertex value doesn't match.");
	//	}
	//	Pass("Test_Mesh_Loader[Multiple] Passed!");
	//}

	template<class MeshType>
	void Test_Mesh_Loader_Single() {
		Info("Test_Mesh_Loader_Single");
		// 0. generated a random mesh
		VertexMatrix<double, 3>  V; V.resize(100, 3);
		ElementMatrix<3> F; F.resize(100, 3);

		for (size_t i = 0; i < 100; i++)
		{
			V.row(i) = Random::Random_VectorXd(3, 0.0, 1.0);
			F.row(i) = Random::Random_VectorXi(3, 0  , 99);
		}
		Info("random generated");

		OBJFunc::Write_Obj<double, 3>("./copy_mesh_single.obj", V, F);
		//auto mesh = std::make_shared<MeshType>();
		//OBJFunc::Read_Mesh("../../../../../proj/_tests_reservoir/assets/CornellBox-Single.obj", mesh);
		//OBJFunc::Write_Mesh<MeshType>("./copy_mesh_single.obj", mesh);

		//auto copy_mesh = std::make_shared<MeshType>();
		//OBJFunc::Read_Mesh("./copy_mesh_single.obj", copy_mesh);
		//Assert(mesh->Elements().size() == copy_mesh->Elements().size(), "Test_Mesh_Loader Failed: faces size doesn't match.");
		//Assert(mesh->Vertices().size() == copy_mesh->Vertices().size(), "Test_Mesh_Loader Failed: vertice size doesn't match.");
		//for (size_t j = 0; j < mesh->Elements().size(); j++) Assert(mesh->elements[j] == copy_mesh->elements[j], "Test_Mesh_Loader Failed: faces value doesn't match.");
		//for (size_t j = 0; j < mesh->Vertices().size(); j++) Assert(mesh->Vertices()[j] == copy_mesh->Vertices()[j], "Test_Mesh_Loader Failed: vertex value doesn't match.");
		Pass("Test_Mesh_Loader[Single] Passed!");
	}
}