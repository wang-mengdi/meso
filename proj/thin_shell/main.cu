#include "AuxFunc.h"
#include "Json.h"
#include "Driver.h"
#include "SoftBodyNonlinearFemThinShell.h"
#include "SoftBodyNonlinearThinShellInitializer.h"
using namespace Meso;

template<int d>
void Run(json& j) {
	SoftBodyNonlinearFemThinShell<d> thin_shell;
	SoftBodyNonlinearThinShellInitializer<d> scene;
	Driver driver;
	driver.Run(j, scene, thin_shell);
}

int main(int argv, char** argc) {
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

		Run<3>(j);
	}
	catch (nlohmann::json::exception& e)
	{
		Meso::Info("json exception {}", e.what());
	}
	return 0;
}