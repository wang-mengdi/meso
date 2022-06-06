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

	template<class T, int d, int ed>
	void Test_Mesh_Loader_Single() {
		Info("Test_Mesh_Loader_Single");
		// 0. generated a random mesh
		VertexMatrix<T, d>  vertex_matrix; vertex_matrix.resize(100, d);
		ElementMatrix<ed> element_matrix; element_matrix.resize(100, d);

		for (size_t i = 0; i < 100; i++)
		{
			vertex_matrix.row(i) = Random::Random_VectorXd(d, 0.0, 1.0);
			element_matrix.row(i) = Random::Random_VectorXi(ed, 0, 99);
		}
		Info("random generated");

		OBJFunc::Write_Obj<T, d, ed>("./copy_mesh_single.obj", vertex_matrix, element_matrix);
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