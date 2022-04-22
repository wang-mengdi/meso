#include "Advection.h"
#include "AuxFunc.h"
#include "Json.h"
#include "FluidEuler.h"
#include "FluidEulerInitializer.h"
#include "Driver.h"
using namespace Meso;

template<int d>
void Run(json &j) {
	FluidEuler<d> fluid;
	FluidEulerInitializer<d> scene;
	Driver driver;
	driver.Run(j, scene, fluid);
}

int main(void) {
	try {
		json j;
		//int dim = 2;
		int dim = Json::Value(j, "dimension", 2);
		if (dim == 2) Run<2>(j);
		else if (dim == 3) Run<3>(j);
	}
	catch (nlohmann::json::exception& e)
	{
		Info("json exception {}", e.what());
	}
	return 0;
}