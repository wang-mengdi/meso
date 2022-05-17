#include "AuxFunc.h"
#include "Json.h"
#include "Driver.h"
#include "SoftBodyNonlinearFemThinShell.h"
#include "SoftBodyNonlinearThinShellInitializer.h"
#include "OptimizerDriver.h"
#include "ThinShellTopologyOptimizer.h"
#include "ThinShellTopologyOptimizerInitializer.h"
#include "omp.h"
using namespace Meso;

template<int d>
void Run(json& j) {
	SoftBodyNonlinearFemThinShell<d> thin_shell;
	SoftBodyNonlinearThinShellInitializer<d> scene;
	Driver driver;
	driver.Run(j, scene, thin_shell);
}

template<int d>
void RunOptimizer(json &j) {
	ThinShellTopologyOptimizer<d> optimizer;
	TopoOptThinShellInitializer<d> scene;
	OptimizerDriver driver;
	driver.Run(j, scene, optimizer);
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

		int thread_num=Json::Value(j, "thread_num", omp_get_max_threads());
		omp_set_num_threads(thread_num);
		
		Info("Using {} threads, out of {} available threads", omp_get_num_threads(), omp_get_max_threads());

		if (j.contains("driver")) {
			Run<3>(j);
		}
		else if (j.contains("optimizer")) {
			RunOptimizer<3>(j);
		}
	}
	catch (nlohmann::json::exception& e)
	{
		Meso::Info("json exception {}", e.what());
	}
	return 0;
}