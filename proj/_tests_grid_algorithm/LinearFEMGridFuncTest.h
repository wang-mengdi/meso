//////////////////////////////////////////////////////////////////////////
// Test linear FEM Grid
// Copyright (c) (2022-), Fan Feng
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "LinearFEMFunc.h"
#include <iostream>
namespace Meso {
	template<int d>
	void Test_Linear_FEM_Grid(real youngs, real poisson) {
			MatrixX K_e;
			MatrixX K_e_table;

			LinearFEMFunc::Cell_Stiffness_Matrix<d>(youngs, poisson, (real)1, K_e);
			LinearFEMFunc::Cell_Stiffness_Matrix<d>(youngs, poisson, K_e_table);
			
			if (!K_e.isApprox(K_e_table)) { 
				Info("K_e:\n{}", K_e);
				Info("K_e_table:\n{}", K_e_table);
				Error("The element stiffness matrix for d={} is not the same as that from table", d);
			}
			else { Pass("The element stiffness matrix for d={} is the same as that from table", d); }
	}
}
