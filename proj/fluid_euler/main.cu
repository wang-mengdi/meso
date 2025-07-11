#include "Advection.h"
#include "AuxFunc.h"
#include "Json.h"
#include "FluidEuler.h"
#include "FluidEulerInitializer.h"
#include "FluidFreeSurface.h"
#include "FluidFreeSurfaceInitializer.h"
#include "Driver.h"
using namespace Meso;

template<int d>
void Run_FluidEuler(json &j) {
	FluidEuler<d> fluid;
	FluidEulerInitializer<d> scene;
	Driver driver;
	driver.Initialize_And_Run(j, scene, fluid);
}


template<int d>
void Run_FluidFreeSurface(json& j) {
	FluidFreeSurface<real, d> fluid;
	FluidFreeSurfaceInitializer<real, d> scene;
	Driver driver;
	driver.Initialize_And_Run(j, scene, fluid);
}

int main(int argc, char **argv) {
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
			Info("Read json file {}", argv[1]);
			std::ifstream json_input(argv[1]);
			json_input >> j;
			json_input.close();
		}

		int dim = Json::Value(j, "dimension", 2);
		std::string simulator=Json::Value(j, "simulator", std::string("euler"));
		
		if (simulator == "euler") {
			if (dim == 2) { Run_FluidEuler<2>(j); }
			else if (dim == 3) { Run_FluidEuler<3>(j); }
		}
		else if (simulator == "freesurface") {
			if (dim == 2) { Run_FluidFreeSurface<2>(j); }
			else if (dim == 3) { Run_FluidFreeSurface<3>(j); }
		}

	}
	catch (nlohmann::json::exception& e)
	{
		Info("json exception {}", e.what());
	}
	return 0;
}