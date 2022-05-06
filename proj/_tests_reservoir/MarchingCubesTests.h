#pragma once

#include "Grid.h"
#include "Field.h"
#include "Mesh.h"
#include "MarchingCubes.h"
namespace Meso {
	void Test_Marching_Cubes() {
		Field<real, 3> f;
		auto m = std::make_shared<TriangleMesh<3>>();
		Marching_Cubes<3>(f, m, 0.);
		return;
	}
}