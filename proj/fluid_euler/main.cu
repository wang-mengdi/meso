#include "Advection.h"
#include "AuxFunc.h"
#include "Json.h"
//#include <nlohmann/json.hpp>
using namespace Meso;


//#include <iostream>
//#include <filesystem>
//namespace fs = std::filesystem;


int main(void) {
	//SemiLagrangian<2> sl;
	real a = MathFunc::Clamp(std::numeric_limits<real>::infinity(), (real)0, (real)1);
	Info("clamped: {}", a);

	//fs::path dir("/tmp");
	return 0;
}