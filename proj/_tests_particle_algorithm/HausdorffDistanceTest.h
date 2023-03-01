#pragma once
#include "HausdorffDistance.h"
#include <iostream>
namespace Meso {
	template<int d,DataHolder side>
	void Test_HausdorffDistance() {
		if constexpr (d == 1) {
			Array<Vector1, side> a(2);
			a[0] = Vector1(1);
			a[1] = Vector1(7);
			Array<Vector1, side> b(2);
			b[0] = Vector1(3);
			b[1] = Vector1(6);
			real result = HausdorffDistance<real,d,side>(a, b);
			Vector1 b_0 = b[0];
			Vector1 a_0 = a[0];
			if (result == (b_0 - a_0).norm()) {
				Pass("Dimension={}, the hausdorff distance test is passed", d);
			}
			else {
				Warn("Dimension={}, the hausdorff distance test is not passed", d);
			}
		}
		else if constexpr (d == 2) {
			Array<Vector2, side> a(4);
			a[0] = Vector2(2, 3);
			a[1] = Vector2(1, 0);
			a[2] = Vector2(4, 8);
			a[3] = Vector2(9, 6);
			Array<Vector2, side> b(3);
			b[0] = Vector2(0, 0);
			b[1] = Vector2(3, 1);
			b[2] = Vector2(5, 4);
			real result = HausdorffDistance<real, d,side>(a,b);
			Vector2 b_2 = b[2];
			Vector2 a_3 = a[3];
			if (result == (a_3-b_2).norm()) {
				Pass("Dimension={}, the hausdorff distance test is passed", d);
			}
			else {
				Warn("Dimension={}, the hausdorff distance test is not passed", d);

			}
		}
		else if constexpr (d == 3) {
			Array<Vector3, side> a(4);
			a[0] = Vector3(2, 3, 0);
			a[1] = Vector3(1, 0, 1);
			a[2] = Vector3(4, 8, 2);
			a[3] = Vector3(9, 6, 3);
			Array<Vector3, side> b(3);
			b[0] = Vector3(0, 0, 2);
			b[1] = Vector3(3, 1, 3);
			b[2] = Vector3(5, 4, 1);

			real result = HausdorffDistance<real, d,side>(a, b);
			Vector3 b_2 = b[2];
			Vector3 a_3 = a[3];
			if (result == (a_3 - b_2).norm()) {
				Pass("Dimension={}, the hausdorff distance test is passed", d);
			}
			else {
				Warn("Dimension={}, the hausdorff distance test is not passed", d);
			}
		}
		else {
			Error("Test_HausdorffDistance: dimension {} not supported",d);
		}
	}
}