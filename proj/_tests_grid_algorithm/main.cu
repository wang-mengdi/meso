#include <iostream>
#include <cstdio>

#include "Common.h"
#include "Grid.h"
#include "Field.h"
#include "LevelsetTest.h"
#include "MarchingCubesTests.h"
#include "LinearFEMGridFuncTest.h"
#include "NonManifoldMarchingCubesTests.h"

using namespace Meso;
#define __THRUST_HAS_CUDART__
int main() {
	//test fast marching
	Test_Fast_Marching<2>(8);
	Test_Fast_Marching<2>(16);
	Test_Fast_Marching<2>(32);
	Test_Fast_Marching<2>(64);
	Test_Fast_Marching<2>(90);
	Test_Fast_Marching<2>(128);

	Test_Fast_Marching<3>(8);
	Test_Fast_Marching<3>(16);
	Test_Fast_Marching<3>(32);
	Test_Fast_Marching<3>(64);
	Test_Fast_Marching<3>(90);
	Test_Fast_Marching<3>(128);

	//test marching cubes
	Test_Marching_Cubes<float>(1, false); // verify by opening marching_cubes.obj
	Test_Marching_Cubes<double>(1, false);

	Test_Non_Manifold_Marching_Cubes1<double>(1, true);
	Test_Non_Manifold_Marching_Cubes2<double>(1, true);
	Test_Non_Manifold_Marching_Cubes3<double>(1, true);

	//test linear FEM Grid
	Test_Linear_FEM_Grid<2>((real)1, 0.3);
	Test_Linear_FEM_Grid<3>((real)1, 0.3);
}