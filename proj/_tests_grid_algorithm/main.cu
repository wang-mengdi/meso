#include <iostream>
#include <cstdio>

#include "Common.h"
#include "Grid.h"
#include "Field.h"
#include "LevelsetPrimeTest.h"
#include "LevelsetAdvanceTest.h"
#include "MarchingCubesTests.h"

using namespace Meso;
#define __THRUST_HAS_CUDART__
int main(){
	//Test_Fast_Marching<3, HOST>(Vector3i(32, 32, 32));
	//return 0;

	Test_Marching_Cubes<float>(1, true); // verify by opening marching_cubes.obj
	Test_Marching_Cubes<double>(1, false);

	Test_Fast_Marching<2, HOST>(Vector2i(192, 168));
	//Test_Fast_Marching<3, HOST>(Vector3i(62, 40, 83));

	Test_Fast_Marching<2, HOST>(Vector2i(8, 8));
	Test_Fast_Marching<2, HOST>(Vector2i(16, 16));
	Test_Fast_Marching<2, HOST>(Vector2i(32, 32));
	Test_Fast_Marching<2, HOST>(Vector2i(64, 64));
	Test_Fast_Marching<2, HOST>(Vector2i(128, 128));

	Test_Fast_Marching<3, HOST>(Vector3i(8, 8, 8));
	Test_Fast_Marching<3, HOST>(Vector3i(16, 16, 16));
	Test_Fast_Marching<3, HOST>(Vector3i(32, 32, 32));
	Test_Fast_Marching<3, HOST>(Vector3i(64, 64, 64));
	Test_Fast_Marching<3, HOST>(Vector3i(128, 128, 128));

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
}