#include <iostream>
#include <cstdio>

#include "Common.h"
#include "Grid.h"
#include "Field.h"
#include "LevelsetPrimeTest.h"
#include "LevelsetAdvanceTest.h"
#include "MarchingCubesTests.h"
#include "LinearFEMGridFuncTest.h"
#include "NonManifoldMarchingCubesTests.h"

using namespace Meso;
#define __THRUST_HAS_CUDART__
int main(){
	//Test_Fast_Marching<2>(8);
	//Test_Fast_Marching<2>(16);
	//Test_Fast_Marching<2>(32);
	//Test_Fast_Marching<2>(64);
	//Test_Fast_Marching<2>(128);

	//Test_Fast_Marching<3>(8);
	//Test_Fast_Marching<3>(16);
	//Test_Fast_Marching<3>(32);
	//Test_Fast_Marching<3>(64);
	//Test_Fast_Marching<3>(128);

	//return 0;

	//Test_Marching_Cubes<float>(1, false); // verify by opening marching_cubes.obj
	//Test_Marching_Cubes<double>(1, false);

	//Test_Non_Manifold_Marching_Cubes1<double>(1, true);
	Test_Non_Manifold_Marching_CubesN<double>(1, true, 64);

	//Test_Fast_Marching<2>(128);
	//Test_Fast_Marching<3>(128);


	// Test_Marching_Cubes<double>(1, false);

	/// prime test
	//Test_Circle_HOST<2, PointIntpLinearPadding0>(0.5, Vector<real,2>::Zero(),  1);
	//Test_Circle_HOST<2, PointIntpLinearPadding0>(0.02, Vector<real, 2>::Zero(), 2);
	//Test_Circle_HOST<2, PointIntpLinearPadding0>(0.5, Vector<real, 2>::Zero()+ Vector<real, 2>(1,2), 1);
	//Test_Circle_HOST<2, PointIntpLinearPadding0>(0.02, Vector<real, 2>::Zero()+Vector<real, 2>(1, 2), 2);


	//Test_Circle_HOST<3, PointIntpLinearPadding0>(0.02, Vector<real, 3>::Zero(), 1);
	//Test_Circle_HOST<3, PointIntpLinearPadding0>(0.2, Vector<real, 3>::Zero(), 2);
	//Test_Circle_HOST<3, PointIntpLinearPadding0>(0.02, Vector<real, 3>::Zero() + Vector<real, 3>(1, 2, 3), 1);
	//Test_Circle_HOST<3, PointIntpLinearPadding0>(0.2, Vector<real, 3>::Zero() + Vector<real, 3>(1, 2, 3), 2);

	//Test_FMM_Circle<2, PointIntpLinearPadding0,HOST>(0.02, Vector<real, 2>::Zero(), 1);
	//Test_FMM_Circle<2, PointIntpLinearPadding0, DEVICE>(0.02, Vector<real, 2>::Zero(), 1);
	//Test_FMM_Circle<3, PointIntpLinearPadding0, HOST>(0.02, Vector<real, 3>::Zero(), 1);
	//Test_FMM_Circle<3, PointIntpLinearPadding0, DEVICE>(0.02, Vector<real, 3>::Zero(), 1);

	//Test_Circle<2, PointIntpLinearPadding0, DEVICE>(0.02, Vector<real, 2>::Zero(), 1);
	//Test_Circle<2, PointIntpLinearPadding0, DEVICE>(0.02, Vector<real, 2>::Zero(), 0.2);
	//Test_Circle<2, PointIntpLinearPadding0, DEVICE>(0.02, Vector<real, 2>::Zero() + Vector<real, 2>(1, 2), 1);
	//Test_Circle<2, PointIntpLinearPadding0, DEVICE>(0.02, Vector<real, 2>::Zero() + Vector<real, 2>(1, 2), 2);


	//Test_Circle<3, PointIntpLinearPadding0, DEVICE>(0.02, Vector<real, 3>::Zero(), 1);
	//Test_Circle<3, PointIntpLinearPadding0, DEVICE>(0.2, Vector<real, 3>::Zero(), 2);
	//Test_Circle<3, PointIntpLinearPadding0, DEVICE>(0.02, Vector<real, 3>::Zero() + Vector<real, 3>(1, 2, 3), 1);
	//Test_Circle<3, PointIntpLinearPadding0, DEVICE>(0.2, Vector<real, 3>::Zero() + Vector<real, 3>(1, 2, 3), 2);


	//Test_Circle<2, PointIntpLinearPadding0, HOST>(0.02, Vector<real, 2>::Zero(), 1);
	//Test_Circle<2, PointIntpLinearPadding0, HOST>(0.02, Vector<real, 2>::Zero(), 0.2);
	//Test_Circle<2, PointIntpLinearPadding0, HOST>(0.02, Vector<real, 2>::Zero() + Vector<real, 2>(1, 2), 1);
	//Test_Circle<2, PointIntpLinearPadding0, HOST>(0.02, Vector<real, 2>::Zero() + Vector<real, 2>(1, 2), 2);


	//Test_Circle<3, PointIntpLinearPadding0, HOST>(0.02, Vector<real, 3>::Zero(), 1);
	//Test_Circle<3, PointIntpLinearPadding0, HOST>(0.2, Vector<real, 3>::Zero(), 2);
	//Test_Circle<3, PointIntpLinearPadding0, HOST>(0.02, Vector<real, 3>::Zero() + Vector<real, 3>(1, 2, 3), 1);
	//Test_Circle<3, PointIntpLinearPadding0, HOST>(0.2, Vector<real, 3>::Zero() + Vector<real, 3>(1, 2, 3), 2);

	//Test_FMM_Circle<2, PointIntpLinearPadding0, HOST>(0.02, Vector<real, 2>::Zero(), 1);
	//Test_FMM_Circle<2, PointIntpLinearPadding0, DEVICE>(0.02, Vector<real, 2>::Zero(), 1);
	//Test_FMM_Circle<3, PointIntpLinearPadding0, HOST>(0.02, Vector<real, 3>::Zero(), 1);
	//Test_FMM_Circle<3, PointIntpLinearPadding0, DEVICE>(0.02, Vector<real, 3>::Zero(), 1);

	//test linear FEM Grid
	//Test_Linear_FEM_Grid<2>((real)1, 0.3);
	//Test_Linear_FEM_Grid<3>((real)1, 0.3);
}