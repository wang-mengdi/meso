#include "NonlinearFemFunc.h"
#include "LinearFemFunc.h"

namespace NonlinearFemFunc
{
	void Set_Dirichlet_Boundary_Helper(SparseMatrix<real>& K, Array<real>& b, const int i, const real psi_D_value)
	{
		for (InnerIterator<real> iter(K, i); iter; ++iter) {
			int j = (int)iter.col();
			if (i == j) { b[j] = psi_D_value * K.coeff(i, j); }
			else {
				real K_ij = K.coeff(i, j); b[j] -= K_ij * psi_D_value;
				K.coeffRef(i, j) = (real)0; K.coeffRef(j, i) = (real)0;
			}
		}
	}

	void D(const Vector<real, 2>& x1, const Vector<real, 2>& x2, Matrix<real, 2>& ds) {
		ds.col(0) = x2 - x1; ds.col(1) = Vector2(x1[1] - x2[1], x2[0] - x1[0]).normalized(); //the perpendicular vector
	}

	void D(const Vector<real, 3>& x1, const Vector<real, 3>& x2, const Vector<real, 3>& x3, Matrix<real, 3>& ds) {
		ds.col(0) = x2 - x1; ds.col(1) = x3 - x1; ds.col(2) = ((x3 - x1).cross(x2 - x1)).normalized(); //the normal vector to the triangle
	}

	void D_Inv_And_Area_And_Normal(const Vector<real, 2>& X1, const Vector<real, 2>& X2,/*rst*/Matrix<real, 2>& dm_inv,/*rst*/real& area, Vector<real, 2>& normal) {
		D(X1, X2, dm_inv); area = (X2 - X1).norm(); normal = dm_inv.col(1); dm_inv = dm_inv.inverse().eval();//area is the segment length
	}

	void D_Inv_And_Area_And_Normal(const Vector<real, 3>& X1, const Vector<real, 3>& X2, const Vector<real, 3>& X3,/*rst*/Matrix<real, 3>& dm_inv,/*rst*/real& area, Vector<real, 3>& normal) {
		D(X1, X2, X3, dm_inv); area = (real)0.5 * ((X2 - X1).cross(X3 - X1)).norm(); normal = dm_inv.col(2); dm_inv = dm_inv.inverse().eval();//area is the triangle area
	}
}