#pragma once
#include "Common.h"
#include "AuxFunc.h"
#include "SparseFunc.h"
#include <iostream>

using namespace Meso;
namespace NonlinearFemFunc
{
	//////////////////////////////////////////////////////////////////////////
	////Thin Shell element

	void D(const Vector<real, 2>& x1, const Vector<real, 2>& x2, Matrix<real, 2>& ds);
	void D(const Vector<real, 3>& x1, const Vector<real, 3>& x2, const Vector<real, 3>& x3, Matrix<real, 3>& ds);

	template<int d> void D(const ArrayF<Vector<real, d>, d>& x, Matrix<real, d>& ds) {
		if constexpr (d == 2) {
			return D(x[0], x[1], ds);
		}
		else if constexpr (d == 3) {
			return D(x[0], x[1], x[2], ds);
		}
		else {
			Assert(false, "dimension {} not supported", d);
		}
	}
	
	template<int d> void D_Inv_And_Area(const Vector<real, d>& X1, const Vector<real, d>& X2, const Vector<real, d>& X3,/*rst*/Matrix<real, d>& dm_inv,/*rst*/real& area) {
		if constexpr (d == 2) {
			D(X1, X2, dm_inv); area = (X2 - X1).norm(); dm_inv = dm_inv.inverse().eval();//area is the triangle area
		}
		else if constexpr (d == 3) {
			D(X1, X2, X3, dm_inv); area = (real)0.5 * ((X2 - X1).cross(X3 - X1)).norm(); dm_inv = dm_inv.inverse().eval();//area is the triangle area
		}
		else {
			Assert(false, "dimension {} not supported", d);
		}
	}

	void D_Inv_And_Area_And_Normal(const Vector<real, 2>& X1, const Vector<real, 2>& X2,/*rst*/Matrix<real, 2>& dm_inv,/*rst*/real& area, Vector<real, 2>& normal);
	void D_Inv_And_Area_And_Normal(const Vector<real, 3>& X1, const Vector<real, 3>& X2, const Vector<real, 3>& X3,/*rst*/Matrix<real, 3>& dm_inv,/*rst*/real& area, Vector<real, 3>& normal);

	template<int d> void D_Inv_And_Area(const ArrayF<Vector<real, d>, d>& x,/*rst*/Matrix<real, d>& dm_inv,/*rst*/real& area) {
		if constexpr (d == 2) {
			D_Inv_And_Area(x[0], x[1], dm_inv, area);
		}
		else if constexpr (d == 3) {
			D_Inv_And_Area(x[0], x[1], x[2], dm_inv, area);
		}
		else {
			Assert(false, "dimension {} not supported", d);
		}
	}

	template<int d> void D_Inv_And_Area_And_Normal(const ArrayF<Vector<real, d>, d>& x,/*rst*/Matrix<real, d>& dm_inv,/*rst*/real& area, Vector<real, d>& normal) {
		if constexpr (d == 2) {
			D_Inv_And_Area_And_Normal(x[0], x[1], dm_inv, area, normal);
		}
		else if constexpr (d == 3) {
			D_Inv_And_Area_And_Normal(x[0], x[1], x[2], dm_inv, area, normal);
		}
		else {
			Assert(false, "dimension {} not supported", d);
		}
	}

	//helper function
	void Set_Dirichlet_Boundary_Helper(SparseMatrix<real>& K, Array<real>& b, const int i, const real psi_D_value);
}