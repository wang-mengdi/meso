#include "Common.h"
#include "HausdorffDistanceTest.h"
using namespace Meso;
int main() {
	Test_HausdorffDistance<1,HOST>();
	Test_HausdorffDistance<2,HOST>();
	Test_HausdorffDistance<3,HOST>();
	Test_HausdorffDistance<1, DEVICE>();
	Test_HausdorffDistance<2, DEVICE>();
	Test_HausdorffDistance<3, DEVICE>();
	return 0;
}