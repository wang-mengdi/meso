#include "SPX_AuxFunc.h"
#include "SPX_Json.h"
#include "SPX_Driver.h"
#include "FluidEulerTwoPhaseClebschDriver.h"

template<int d>
void Run(json& j) {
	FluidEulerTwoPhaseClebschDriver<d> driver;
	driver.scale = Json::Value<int>(j, "scale", 64);
	driver.output_dir = Json::Value<std::string>(j, "output_dir", "out");
	driver.test = Json::Value<int>(j, "test", 1);
	driver.last_frame = Json::Value<int>(j, "last_frame", 500);
	driver.frame_rate = Json::Value<int>(j, "frame_rate", 20);
	driver.cfl = Json::Value<real>(j, "cfl", 1.);
	driver.beta = Json::Value<real>(j, "beta", 0.95);
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
		//Meso::Info("json exception {}", e.what());
		fmt::print("json exception {}", e.what());
	}
	return 0;
}