#include "SPX_AuxFunc.h"
#include "SPX_Json.h"
#include "SPX_Driver.h"
#include "FluidEulerTwoPhaseClebschDriver.h"

template<int d>
void Run(json& j) {
	j = j.at("driver");
	FluidEulerTwoPhaseClebschDriver<d> driver;
	driver.scale = Json::Value<int>(j, "scale", 64);
	driver.output_dir = Json::Value<std::string>(j, "output_dir", "output");
	driver.test = Json::Value<int>(j, "test", 1);
	driver.last_frame = Json::Value<int>(j, "last_frame", 1000);
	driver.frame_rate = Json::Value<int>(j, "fps", 20);
	driver.cfl = Json::Value<real>(j, "cfl", 1.);
	driver.rho_A = Json::Value<real>(j, "rho_A", 0.001);
	driver.rho_L = Json::Value<real>(j, "rho_L", 1.);
	driver.h_bar = Json::Value<real>(j, "h_bar", 0.01);
	driver.beta = Json::Value<real>(j, "beta", 0.95);
	driver.Initialize();
	driver.Run();
}

int main(int argv, char** argc) {
	//Info("omp max num threads: {}", omp_get_max_threads());

	try {
		json j = {
			{
				"driver",
				{
					{"last_frame",10}
				}
			},
		};
		if (argv > 1) {
			std::ifstream json_input(argc[1]);
			json_input >> j;
			json_input.close();
		}
		int dim = Json::Value(j, "dimension", 2);
		if (dim == 2) Run<2>(j);
		else if (dim == 3) Run<3>(j);
	}
	catch (nlohmann::json::exception& e)
	{
		fmt::print("json exception {}", e.what());
	}
	return 0;
}