#include "Json.h"
#include "FluidPIC.h"
#include "FluidPICInitializer.h"
#include "Driver.h"
using namespace Meso;

template<int d>
void Run_FluidPIC(json& j) {
	FluidPIC<d> fluid;
	FluidPICInitializer<d> scene;
	Driver driver;
	driver.Run(j, scene, fluid);
}

int main(int argc, char** argv) {
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
		if (argc > 1) {
			std::ifstream json_input(argv[1]);
			json_input >> j;
			json_input.close();
		}

		int dim = Json::Value(j, "dimension", 2);
		std::string simulator = Json::Value(j, "simulator", std::string("pic"));

		if (simulator == "pic") {
			if (dim == 2) { Run_FluidPIC<2>(j); }
			else if (dim == 3) { Run_FluidPIC<3>(j); }
		}
	}
	catch (nlohmann::json::exception& e)
	{
		Info("json exception {}", e.what());
	}
	return 0;
}