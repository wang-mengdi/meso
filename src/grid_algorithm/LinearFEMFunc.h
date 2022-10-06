#pragma once
#include "Common.h"
#include "Constants.h"
#include "SparseFunc.h"


namespace Meso {
	namespace LinearFEMFunc
	{
		const static Vector3i CORNER_OFFSETS[] = {Vector3i(0,0,0), Vector3i(1,0,0), Vector3i(1,1,0), Vector3i(0,1,0), Vector3i(0,0,1), Vector3i(1,0,1), Vector3i(1,1,1), Vector3i(0,1,1) };
		//correspondence of index to corner coordinate given a center coordinate
		template<int d> Vector<int, d> Corner_Offset(const Vector<int, d>& center, int i) {
			Assert(i < pow(2, d), "Corner_Offset: index out of  range");
			if constexpr (d == 2) { return center + Vector2i(i & 0x1, (i >> 1) & 0x1); }
			else if constexpr (d == 3) { return center + CORNER_OFFSETS[i];}
			else { Error("Corner_Offset: dimension not supported"); return Vector<int, d>(); }
		}

		////////////////////////////////////////////////////////////////////////
		// E: strain -> stress
		template<int d> void Strain_Stress_Matrix_Linear(const real youngs, const real poisson, MatrixX& E);

		////////////////////////////////////////////////////////////////////////
		//Hex element
		template<int d> void dNde(const Vector<real, d>& natural_coord, MatrixX& dNde);							//dNde, return a dxn matrix, n is the number of basis functions
		template<int d> void Cell_dNdX(const Vector<real, d>& natural_coord, const real dx, MatrixX& dNdX);		//dNdX

		template<int d> void Cell_Strain_Displacement_Matrix(const Vector<real, d>& natural_coord, const real dx, MatrixX& B);			//B: displacement->strain
		template<int d> void Cell_Stiffness_Matrix(const real youngs, const real poisson, const real dx, MatrixX& K_e);
		template<int d> void Cell_Stiffness_Matrix(const real youngs, const real poisson, MatrixX& K_e);								//Use a table to caculate the stiffness matrix

		//Gaussian integration
		template<int d> void Initialize_Gaussian_Integration_Points(Matrix<real,1<<d,d>& points, Vector<real, 1<<d>& weights); // Gauss product rule p = 2, 3rd order accuracy

		////////////////////////////////////////////////////////////////////////
		//Operations on the global stiffness matrix
		template<int d> void Add_Cell_Stiffness_Matrix(SparseMatrix<real>& K, const MatrixX& K_e, const Array<int>& nodes);

		void Set_Dirichlet_Boundary_Helper(SparseMatrix<real>& K, VectorX& b, const int i, const real psi_D_value);
	}
}