#include "AuxFunc.h"
#include "Json.h"
#include "omp.h"
#include "Driver.h"
#include "DiscreteShell.h"
#include "DiscreteShellInitializer.h"
#include "OptimizerDriver.h"
#include "DiscreteShellQuasistatic.h"
#include "DiscreteShellTopologyOptimizer.h"
#include "DiscreteShellTopoOptInitializer.h"
using namespace Meso;

template<int d>
void Run(json& j) {
	DiscreteShell<d> thin_shell;
	DiscreteShellInitializer<d> scene;
	Driver driver;
	driver.Run(j, scene, thin_shell);
}

template<int d>
void RunQuasistaticSolver(json& j) {
	DiscreteShellQuasistatic<d> optimizer;
	DiscreteShellInitializer<d> scene;
	OptimizerDriver driver;
	driver.Run(j, scene, optimizer);
}

//template<int d>
//void RunOptimizer(json &j) {
//	DiscreteShellTopologyOptimizer<d> optimizer;
//	DiscreteShellTopoOptInitializer<d> scene;
//	OptimizerDriver driver;
//	driver.Run(j, scene, optimizer);
//}

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
			RunQuasistaticSolver<3>(j);
		}
	}
	catch (nlohmann::json::exception& e)
	{
		Meso::Info("json exception {}", e.what());
	}
	return 0;
}