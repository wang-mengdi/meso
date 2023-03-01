#include "Common.h"
#include "HausdorffDistanceTest.h"
using namespace Meso;
int main() {
	Test_HausdorffDistance<1>();
	Test_HausdorffDistance<2>();
	Test_HausdorffDistance<3>();
	return 0;
}