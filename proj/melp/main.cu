#include "AuxFunc.h"
#include "Json.h"
#include "IOFunc.h"
#include "Driver.h"
#include "MELP.h"
#include "MELPInitializer.h"
#include "Points.h"
#include "NAParticles.h"
#include "Common.h"
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
	NAParticles<d> naps;
	//naps.Add_Attribute("s", (float)10);
	//naps.Add_Attribute("t", VectorFunc::V<d>(1,1,1));
	naps.template Add_Attribute<float>("s", 10);
	naps.template Add_Attribute<Vector3>("t", Vector3::Ones());
	naps.Resize(2);
	Info("size: {}", naps.Size());
	//legal
	Array<float>& data = naps.template Get_Attribute<float>("s");
	Info("size 2: {}", data.size());
	Info("data 0: {}", naps.template Get_Entry<float>("s", 0));
	Info("data 1: {}", naps.template Get_Entry<float>("s", 1));
	Info("X 0: {}", naps.Get["x"](0));
	Info("X 1: {}", naps.Get["x"](1));
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