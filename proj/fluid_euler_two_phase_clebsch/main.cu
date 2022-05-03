#include "SPX_AuxFunc.h"
#include "Json.h"
#include "SPX_Driver.h"
#include "FluidEulerTwoPhaseClebschDriver.h"

template<int d>
void Run(json& j) {
	FluidEulerTwoPhaseClebschDriver<d> driver;
	driver.scale = j.get<int>("scale", 64);
	driver.output_dir = j.get<std::string>("output_dir", "out");
	driver.test = j.get<int>("test", 1);
	driver.last_frame = j.get<int>("last_frame", 500);
	driver.frame_rate = j.get<int>("frame_rate", 20);
	driver.cfl = j.get<real>("cfl", 1.);
	driver.beta = j.get<real>("beta", 0.95);
	driver.Initialize();
	driver.Run();
}

int main(int argv, char** argc) {
	try {
		json j;
		if (argv > 1) {
			std::ifstream json_input(argc[1]);
			json_input >> j;
			json_input.close();
		}

		Run<3>(j);
	}
	catch (nlohmann::json::exception& e)
	{
		Meso::Info("json exception {}", e.what());
	}
	return 0;
}