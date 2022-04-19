#include "Advection.h"
#include "AuxFunc.h"
#include "Json.h"
//#include <nlohmann/json.hpp>
using namespace Meso;


//#include <iostream>
//#include <filesystem>
//namespace fs = std::filesystem;

void foo(const std::string s) {
	Info("foo {}", s);
}

int main(void) {
	//SemiLagrangian<2> sl;
	real a = MathFunc::Clamp(std::numeric_limits<real>::infinity(), (real)0, (real)1);
	std::string c = "123";
	Info("clamped: {} {}", a, c);
	foo("345");

	Check_Cuda_Memory("main");

	//fs::path dir("/tmp");
	return 0;
}