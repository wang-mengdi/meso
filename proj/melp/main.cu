#include "AuxFunc.h"
#include "Json.h"
#include "IOFunc.h"
#include "Driver.h"
#include "MELP.h"
#include "MELPInitializer.h"
#include "Points.h"
#include "NAParticles.h"
#include "EulerParticles.h"
#include "Common.h"
#include "PointsCreator.h"

using namespace Meso;

template<int d>
void Run(json &j) {
	std::cout << "Running" << std::endl;
	MELP<d> fluid;
	MELPInitializer<d> scene;
	Driver driver;
	driver.Initialize_And_Run(j, scene, fluid);
}

template<int d>
void Test_Points(void) {
}

template<int d>
void Test_NAParticles(void) {
	Typedef_VectorD(d);
	NAParticles<d> naps;
	if constexpr (d == 3) {
		Initialize_Lattice_Points(Vector3::Zero(), 3, 3, Vector3::Unit(0), Vector3::Unit(1), naps, naps.xRef());
		//Initialize_Sphere_Points_Regular(Vector3::Zero(), 1, 100, naps, naps.xRef());
	}
	naps.Update_Searcher();

	Array<int> nbs;
	naps.nbs_searcher->Find_Nbs(Vector3::Zero(), 1.5, nbs);
	Info("num nbs: {}",nbs.size());
	naps.nbs_searcher->Find_Nbs(Vector3::Zero(), 15, nbs);
	Info("num nbs 2: {}", nbs.size());
}

template<int d>
void Test_EParticles(void) {
	Typedef_VectorD(d);
	EulerParticles<d> eps;
	//if constexpr (d == 3) {
	//	Initialize_Lattice_Points(Vector3::Zero(), 3, 3, 1, 1, eps, "x");
	//}
	//eps.Update_Searcher();
	//Array<int> nbs;
	//eps.nbs_searcher->Find_Neighbors(Vector3::Zero(), 1.5, nbs);
	//Info("num nbs: {}", nbs.size());
	//eps.nbs_searcher->Find_Neighbors(Vector3::Zero(), 15, nbs);
	//Info("num nbs 2: {}", nbs.size());
	//Info("E0: {}", eps.E(0).col(0));
}

template<int d>
void Test_MELP(void) {
	Typedef_VectorD(d);
}

template<int d>
void Test_CG_Solver(void) {
	//Typedef_VectorD(d);
	//Test_Sparse_Matrix();
}

int main(int argv, char **argc) {
	try {
		json j = {
			{
				"driver",
				{
					{"last_frame",10}
				}
			},
			{"scene",json::object()}
		};
		if (argv > 1) {
			std::ifstream json_input(argc[1]);
			json_input >> j;
			json_input.close();
		}
		int dim = Json::Value(j, "dimension", 2);
		if (dim == 2) Run<2>(j);
		else if (dim == 3) Run<3>(j);
		//Test_NAParticles<3>();
		//Test_EParticles<3>();
		//Test_CG_Solver<3>();
	}
	catch (nlohmann::json::exception& e)
	{
		Info("json exception {}", e.what());
	}
	return 0;
}