#include <iostream>
#include <cstdio>

#include "Common.h"
#include "Grid.h"
#include "Field.h"
#include "LevelsetPrimeTest.h"
#include "LevelsetAdvanceTest.h"

using namespace Meso;
#define __THRUST_HAS_CUDART__
int main(){
	Test_Marching_Cubes<float, 3>(); // verify by opening marching_cubes.obj

	/// prime test
	//Test_Circle_HOST<2, PointIntpLinearPadding0>(0.5, Vector<real,2>::Zero(),  1);
	//Test_Circle_HOST<2, PointIntpLinearPadding0>(0.02, Vector<real, 2>::Zero(), 2);
	//Test_Circle_HOST<2, PointIntpLinearPadding0>(0.5, Vector<real, 2>::Zero()+ Vector<real, 2>(1,2), 1);
	//Test_Circle_HOST<2, PointIntpLinearPadding0>(0.02, Vector<real, 2>::Zero()+Vector<real, 2>(1, 2), 2);


	//Test_Circle_HOST<3, PointIntpLinearPadding0>(0.02, Vector<real, 3>::Zero(), 1);
	//Test_Circle_HOST<3, PointIntpLinearPadding0>(0.2, Vector<real, 3>::Zero(), 2);
	//Test_Circle_HOST<3, PointIntpLinearPadding0>(0.02, Vector<real, 3>::Zero() + Vector<real, 3>(1, 2, 3), 1);
	//Test_Circle_HOST<3, PointIntpLinearPadding0>(0.2, Vector<real, 3>::Zero() + Vector<real, 3>(1, 2, 3), 2);

	Test_FMM_Circle<2, PointIntpLinearPadding0,HOST>(0.02, Vector<real, 2>::Zero(), 1);
	//Test_FMM_Circle<2, PointIntpLinearPadding0, DEVICE>(0.02, Vector<real, 2>::Zero(), 1);
	//Test_FMM_Circle<3, PointIntpLinearPadding0, HOST>(0.02, Vector<real, 3>::Zero(), 1);
	//Test_FMM_Circle<3, PointIntpLinearPadding0, DEVICE>(0.02, Vector<real, 3>::Zero(), 1);

	Test_Circle<2, PointIntpLinearPadding0, DEVICE>(0.02, Vector<real, 2>::Zero(), 1);
	Test_Circle<2, PointIntpLinearPadding0, DEVICE>(0.02, Vector<real, 2>::Zero(), 0.2);
	Test_Circle<2, PointIntpLinearPadding0, DEVICE>(0.02, Vector<real, 2>::Zero() + Vector<real, 2>(1, 2), 1);
	Test_Circle<2, PointIntpLinearPadding0, DEVICE>(0.02, Vector<real, 2>::Zero() + Vector<real, 2>(1, 2), 2);


	Test_Circle<3, PointIntpLinearPadding0, DEVICE>(0.02, Vector<real, 3>::Zero(), 1);
	Test_Circle<3, PointIntpLinearPadding0, DEVICE>(0.2, Vector<real, 3>::Zero(), 2);
	Test_Circle<3, PointIntpLinearPadding0, DEVICE>(0.02, Vector<real, 3>::Zero() + Vector<real, 3>(1, 2, 3), 1);
	Test_Circle<3, PointIntpLinearPadding0, DEVICE>(0.2, Vector<real, 3>::Zero() + Vector<real, 3>(1, 2, 3), 2);


	Test_Circle<2, PointIntpLinearPadding0, HOST>(0.02, Vector<real, 2>::Zero(), 1);
	Test_Circle<2, PointIntpLinearPadding0, HOST>(0.02, Vector<real, 2>::Zero(), 0.2);
	Test_Circle<2, PointIntpLinearPadding0, HOST>(0.02, Vector<real, 2>::Zero() + Vector<real, 2>(1, 2), 1);
	Test_Circle<2, PointIntpLinearPadding0, HOST>(0.02, Vector<real, 2>::Zero() + Vector<real, 2>(1, 2), 2);


	Test_Circle<3, PointIntpLinearPadding0, HOST>(0.02, Vector<real, 3>::Zero(), 1);
	Test_Circle<3, PointIntpLinearPadding0, HOST>(0.2, Vector<real, 3>::Zero(), 2);
	Test_Circle<3, PointIntpLinearPadding0, HOST>(0.02, Vector<real, 3>::Zero() + Vector<real, 3>(1, 2, 3), 1);
	Test_Circle<3, PointIntpLinearPadding0, HOST>(0.2, Vector<real, 3>::Zero() + Vector<real, 3>(1, 2, 3), 2);

	Test_FMM_Circle<2, PointIntpLinearPadding0,HOST>(0.02, Vector<real, 2>::Zero(), 1);
	Test_FMM_Circle<2, PointIntpLinearPadding0, DEVICE>(0.02, Vector<real, 2>::Zero(), 1);
	Test_FMM_Circle<3, PointIntpLinearPadding0, HOST>(0.02, Vector<real, 3>::Zero(), 1);
	Test_FMM_Circle<3, PointIntpLinearPadding0, DEVICE>(0.02, Vector<real, 3>::Zero(), 1);
}