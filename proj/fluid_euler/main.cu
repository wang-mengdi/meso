#include "Advection.h"
#include "AuxFunc.h"
#include "Json.h"
#include "FluidEuler.h"
#include "FluidEulerInitializer.h"
#include "FluidImpulse.h"
#include "FluidImpulseInitializer.h"
#include "FluidFreeSurface.h"
#include "FluidFreeSurfaceInitializer.h"
#include "Driver.h"
using namespace Meso;

template<int d>
void Run_Fluid_Euler(json &j) {
	FluidEuler<d> fluid;
	FluidEulerInitializer<d> scene;
	Driver driver;
	driver.Run(j, scene, fluid);
}

template<int d>
void Run_Fluid_Impulse(json& j) {
	FluidImpulse<d> fluid;
	FluidImpulseInitializer<d> scene;
	Driver driver;
	driver.Run(j, scene, fluid);
}

int main(int argc, char **argv) {
	FluidFreeSurface<real, 2> fsf;
	FluidFreeSurfaceInitializer<real, 2> fit;
	return 0;

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
		std::string simulator=Json::Value(j, "simulator", std::string("euler"));
		
		if (simulator == "euler") {
			if (dim == 2) { Run_Fluid_Euler<2>(j); }
			else if (dim == 3) { Run_Fluid_Euler<3>(j); }
		}
		else if (simulator == "impulse") {
			if (dim == 2) { Run_Fluid_Impulse<2>(j); }
			else if (dim == 3) { Run_Fluid_Impulse<3>(j); }
		}

	}
	catch (nlohmann::json::exception& e)
	{
		Info("json exception {}", e.what());
	}
	return 0;
}