#include "AuxFunc.h"
#include "Json.h"
#include "IOFunc.h"
#include "Driver.h"
#include "MELP.h"
#include "MELPInitializer.h"
#include "Points.h"
#include "NAParticles.h"
#include "Common.h"
#include "InitializePoints.h"

using namespace Meso;

template<int d>
void Run(json &j) {
	std::cout << "Running" << std::endl;
	MELP<d> fluid;
	MELPInitializer<d> scene;
	Driver driver;
	driver.Run(j, scene, fluid);
}

template<int d>
void Test_Points(void) {
	//Typedef_VectorD(d)
	//std::cout << "Testing Points" << std::endl;
	//Points ps;
	//ps.Add_Attribute<Vector3>("x");
	//ps.Add_Attribute<Vector2>("v");
	//ps.Resize(123);
	//std::cout << "map size: " << ps.att_map.size() << std::endl;
	//std::cout << "size: " << ps.Size() << std::endl;
	//std::cout << "Done Testing Points" << std::endl;
}

template<int d>
void Test_NAParticles(void) {
	Typedef_VectorD(d);
	NAParticles<d> naps;
	naps.Resize(2);
	Info("shit: {}", naps.x(0));
	Info("data 0: {}", naps.template Get_Entry<VectorD>("x", 0));
	Array<VectorD> ref = naps.xRef();
	Info("shit1: {}", ref[0]);
	Info("shit2: {}", ref[1]);
	//naps.Resize(2);
	//Info("size: {}", naps.Size());
	////legal
	//Array<float>& data = naps.template Get_Attribute<float>("s");
	//Info("size 2: {}", data.size());
	//Info("data 0: {}", naps.template Get_Entry<float>("s", 0));
	//Info("data 1: {}", naps.template Get_Entry<float>("s", 1));
	//Info("X 0: {}", naps.Get["x"](0));
	//Info("X 1: {}", naps.Get["x"](1));

	//if constexpr (d == 3) {
	//	Initialize_Lattice_Points(Vector3::Zero(), 3, 3, 1, 1, naps, "x");
	//}
	////Array<VectorD> positions = naps.Ref["x"]();
	////for (int i = 0; i < positions.size(); i++) {
	////	std::cout << "Pos: \n" << positions[i] << std::endl;
	////}
	//naps.Update_Searcher();
	//Array<int> nbs;
	//naps.nbs_searcher->Find_Neighbors(Vector3::Zero(), 1.5, nbs);
	//Info("num nbs: {}",nbs.size());
	//naps.nbs_searcher->Find_Neighbors(Vector3::Zero(), 15, nbs);
	//Info("num nbs 2: {}", nbs.size());
	//Info("data 2: {}", naps.template Get<float>("s", 2)); //illegal
	//Info("shitman: {}", shitman(0.5,3.7));
	//illegal
	//Array<Vector3>& data2 = naps.template Get_Attribute<Vector3>("s");
	//Info("size 2: {}", data2.size());
	//Info("data 1: {}", data2[0]);
	//Info("data 2: {}", data2[1]);
}

//template<int d>
//void Test_NAParticles_2(void) {
//	NAParticles<d> naps;
//	Register_Attribute(naps, real, "s", 10);
//	naps.Resize(2);
//	Info("size: {}", naps.Size());
//	//legal
//	Array<real>& data = naps.template Get_Attribute<real>("s");
//	Info("size 2: {}", data.size());
//	Info("data 0: {}", naps.template Get<real>("s", 0));
//	Info("data 1: {}", naps.template Get<real>("s", 1));
//}

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
		//int dim = 2;
		//int dim = Json::Value(j, "dimension", 2);
		//if (dim == 2) Run<2>(j);
		//else if (dim == 3) Run<3>(j);
		//Test_Points<3>();
		//Test_NAParticles();
		Test_NAParticles<3>();
	}
	catch (nlohmann::json::exception& e)
	{
		Info("json exception {}", e.what());
	}
	return 0;
}