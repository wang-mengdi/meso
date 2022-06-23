#include "Grid.h"
#include "AuxFunc.h"
#include "LinearFEMFunc.h"

namespace Meso {
	namespace LinearFEMFunc {
		////////////////////////////////////////////////////////////////////////
		//Material model
		template<> void Strain_Stress_Matrix_Linear<2>(const real youngs, const real poisson, MatrixX& E)
		{
			////zero stress in z
			E.resize(3, 3); E.fill((real)0); real e = youngs / ((real)1 - pow(poisson, 2)); real G = youngs / ((real)2 * ((real)1 + poisson));
			E(0, 0) = E(1, 1) = e; E(0, 1) = E(1, 0) = e * poisson; E(2, 2) = G;
		}

		template<> void Strain_Stress_Matrix_Linear<3>(const real youngs, const real poisson, MatrixX& E)
		{
			E.resize(6, 6); E.fill((real)0); real e = youngs / (((real)1 - (real)2 * poisson) * ((real)1 + poisson)); real G = youngs / ((real)2 * ((real)1 + poisson));
			E(0, 0) = E(1, 1) = E(2, 2) = e * ((real)1 - poisson); E(0, 1) = E(1, 0) = E(0, 2) = E(2, 0) = E(1, 2) = E(2, 1) = e * poisson; E(3, 3) = E(4, 4) = E(5, 5) = G;
		}

		////////////////////////////////////////////////////////////////////////
		//Hex element
		template<> void dNde<2>(const Vector2& natural_coord, MatrixX& dNde)
		{	//N0=1/4(1-s)(1-t)
			//N1=1/4(1+s)(1-t)
			//N2=1/4(1-s)(1+t)
			//N3=1/4(1+s)(1+t)
			
			/*
			dN1/ds dN2/ds dN3/ds dN3/ds
			dN1/dt dN2/dt dN3/dt dN4/dt
			*/

			/* the order of the shape function matters
			* how much each node infuence the result
			N2	N3
			N0	N1
			*/
			
			const real s = natural_coord[0]; const real t = natural_coord[1];
			dNde << t - 1,  1 - t,  -1 - t, t + 1,
					s - 1, -1 - s,   1 - s, s + 1;
			dNde *= (real).25;
		}

		template<> void dNde<3>(const Vector3& natural_coord, MatrixX& dNde)
		{
			const real s = natural_coord[0]; const real t = natural_coord[1]; const real q = natural_coord[2];
			dNde << -(1 - t) * (1 - q), -(1 - t) * (1 + q), -(1 + t) * (1 - q), -(1 + t) * (1 + q), (1 - t)* (1 - q), (1 - t)* (1 + q), (1 + t)* (1 - q), (1 + t)* (1 + q),
					-(1 - s) * (1 - q), -(1 - s) * (1 + q), (1 - s)* (1 - q), (1 - s)* (1 + q), -(1 + s) * (1 - q), -(1 + s) * (1 + q), (1 + s)* (1 - q), (1 + s)* (1 + q),
					-(1 - s) * (1 - t), (1 - s)* (1 - t), -(1 - s) * (1 + t), (1 - s)* (1 + t), -(1 + s) * (1 - t), (1 + s)* (1 - t), -(1 + s) * (1 + t), (1 + s)* (1 + t);
			dNde *= (real).125;
		}

		template<int d> void Cell_dNdX(const Vector<real, d>& natural_coord, const real dx, MatrixX& dNdX)
		{
			dNde<d>(natural_coord, dNdX); dNdX *= ((real)2 / dx); //changed from 0.5/dx
		}

		template<> void Cell_Strain_Displacement_Matrix<2>(const Vector2& natural_coord, const real dx, MatrixX& B)
		{	
			/*
				dN0/dx		0	dN1/dx	0		dN2/dx		0	dN3/dx	0
				0		dN0/dy	0		dN1/dy	0		dN2/dy	0		dN3/dy	
				dN0/dy	dN0/dx	dN1/dy	dN1/dx	N2/dy	dN2/dx	dN3/dy	dN3/dx
			*/

			int vtx_n = (int)pow(2, 2);
			B.resize(3, 8); B.fill(0);
			MatrixX dNdX(2, vtx_n); Cell_dNdX<2>(natural_coord, dx, dNdX);
			const int x = 0; const int y = 1;
			B(0, 0) = dNdX(x, 0); B(0, 2) = dNdX(x, 1); B(0, 4) = dNdX(x, 2); B(0, 6) = dNdX(x, 3);
			B(1, 1) = dNdX(y, 0); B(1, 3) = dNdX(y, 1);	B(1, 5) = dNdX(y, 2); B(1, 7) = dNdX(y, 3);
			B(2, 0) = dNdX(y, 0); B(2, 1) = dNdX(x, 0);	B(2, 2) = dNdX(y, 1); B(2, 3) = dNdX(x, 1);	
			B(2, 4) = dNdX(y, 2); B(2, 5) = dNdX(x, 2);	B(2, 6) = dNdX(y, 3); B(2, 7) = dNdX(x, 3);
		}

		template<> void Cell_Strain_Displacement_Matrix<3>(const Vector3& natural_coord, const real dx, MatrixX& B)
		{
			const int d = 3; int tensor_n = d * (d + 1) / 2; int r = tensor_n; int vtx_n = (int)pow(2, d); int c = d * vtx_n;
			B.resize(r, c); B.fill(0);
			MatrixX dNdX(d, vtx_n); Cell_dNdX<d>(natural_coord, dx, dNdX);

			const int x = 0; const int y = 1; const int z = 2;
			Set_Cell_B_Elements_Helper<d>(0, 0, dNdX.row(x), B);
			Set_Cell_B_Elements_Helper<d>(1, 1, dNdX.row(y), B);
			Set_Cell_B_Elements_Helper<d>(2, 2, dNdX.row(z), B);
			Set_Cell_B_Elements_Helper<d>(3, 0, dNdX.row(y), B);
			Set_Cell_B_Elements_Helper<d>(3, 1, dNdX.row(x), B);
			Set_Cell_B_Elements_Helper<d>(4, 1, dNdX.row(z), B);
			Set_Cell_B_Elements_Helper<d>(4, 2, dNdX.row(y), B);
			Set_Cell_B_Elements_Helper<d>(5, 0, dNdX.row(z), B);
			Set_Cell_B_Elements_Helper<d>(5, 2, dNdX.row(x), B);
		}

		template<int d> void Cell_Stiffness_Matrix(const real youngs, const real poisson, const real dx, MatrixX& K_e)
		{
			int n = d * (int)pow(2, d); K_e.resize(n, n); K_e.fill((real)0);
			ArrayF2P<Vector<real, d>, d> points; ArrayF2P<real, d> weights; 
			Initialize_Gaussian_Integration_Points<d>(points, weights);
			MatrixX E; Strain_Stress_Matrix_Linear<d>(youngs, poisson, E);
			real J_det = pow(dx / (real)2, (real)d);

			for (auto i = 0; i < points.size(); i++) {
				MatrixX B; Cell_Strain_Displacement_Matrix<d>(points[i], dx, B);	
				MatrixX K0; K0 = B.transpose() * E * B * J_det * weights[i];
				K_e += K0;
			}
		}
		template void Cell_Stiffness_Matrix<2>(const real youngs, const real poisson, const real dx, MatrixX& K_e);
		template void Cell_Stiffness_Matrix<3>(const real youngs, const real poisson, const real dx, MatrixX& K_e);

		template<> void Cell_Stiffness_Matrix<2>(const real youngs, const real poisson, MatrixX& K_e)
		{
			int n = 8; K_e.resize(n, n); K_e.fill((real)0);
			/*1 / 2 - nu / 6, 1 / 8 + nu / 8, - 1 / 4 - nu / 12, - 1 / 8 + 3 * nu / 8
			- 1 / 4 + nu / 12, - 1 / 8 - nu / 8, nu / 6, 1 / 8 - 3 * nu / 8*/
			VectorX k(n);
			k << (real)1 / (real)2 - poisson / (real)6, (real)1 / (real)8 + poisson / (real)8, -(real)1 / (real)4 - poisson / (real)12, -(real)1 / (real)8 + (real)3 * poisson / (real)8,
				- (real)1 / (real)4 + poisson / (real)12, -(real)1 / (real)8 - poisson / (real)8, poisson / (real)6, (real)1 / (real)8 - (real)3 * poisson / (real)8;
			//This is different than the table in 99 topo because we have different order of shape function: N2 and N3 are switched
			/*[k(1) k(2) k(3) k(4) k(7) k(8) k(5) k(6) 
			   k(2) k(1) k(8) k(7) k(4) k(3) k(6) k(5) 
			   k(3) k(8) k(1) k(6) k(5) k(2) k(7) k(4) 
			   k(4) k(7) k(6) k(1) k(2) k(5) k(8) k(3) 
			   k(7) k(4) k(5) k(2) k(1) -k(2) k(3) -k(4)
			   k(8) k(3) k(2) k(5) -k(2) k(1) -k(8) k(7)
			   k(5) k(6) k(7) k(8) k(3) -k(8) k(1) -k(6)
			   k(6) k(5) k(4) k(3) -k(4) k(7) -k(6) k(1)]*/

			K_e.row(0) << k[0], k[1], k[2], k[3], k[6], k[7], k[4], k[5];
			K_e.row(1) << k[1], k[0], k[7], k[6], k[3], k[2], k[5], k[4];
			K_e.row(2) << k[2], k[7], k[0], k[5], k[4], k[1], k[6], k[3];
			K_e.row(3) << k[3], k[6], k[5], k[0], k[1], k[4], k[7], k[2];
			K_e.row(4) << k[6], k[3], k[4], k[1], k[0], -k[1], k[2], -k[3];
			K_e.row(5) << k[7], k[2], k[1], k[4], -k[1], k[0], -k[7], k[6];
			K_e.row(6) << k[4], k[5], k[6], k[7], k[2], -k[7],  k[0], -k[5];
			K_e.row(7) << k[5], k[4], k[3], k[2], -k[3], k[6], -k[5], k[0];

			K_e *= youngs;
			K_e /= ((real)1 - pow(poisson, (real)2));
		}

		template<> void Cell_Stiffness_Matrix<3>(const real youngs, const real poisson, MatrixX& K_e)
		{
			Info("Not implemented yet");
			return;
		}

		template<int d> void Set_Cell_B_Elements_Helper(const int r, const int c, const VectorX& dN, MatrixX& B)
		{
			for (int i = 0; i < (int)dN.size(); i++)B(r, c + i * d) = dN[i];
		}

		////////////////////////////////////////////////////////////////////////
		//Gaussian integration
		template<> void Initialize_Gaussian_Integration_Points<2>(ArrayF2P<Vector2, 2>& points, ArrayF2P<real, 2>& weights) //2^2=4 sample points
		{
			real c = (real)1 / sqrt((real)3); points[0] = Vector2(-c, -c); points[1] = Vector2(-c, c); points[2] = Vector2(c, -c); points[3] = Vector2(c, c); ArrayFunc::Fill(weights, (real)1);
		}

		template<> void Initialize_Gaussian_Integration_Points<3>(ArrayF2P<Vector3, 3>& points, ArrayF2P<real, 3>& weights) //2^3=8 sample points
		{
			real c = (real)1 / sqrt((real)3); points[0] = Vector3(-c, -c, -c); points[1] = Vector3(-c, -c, c); points[2] = Vector3(-c, c, -c); points[3] = Vector3(-c, c, c);
			points[4] = Vector3(c, -c, -c); points[5] = Vector3(c, -c, c); points[6] = Vector3(c, c, -c); points[7] = Vector3(c, c, c); ArrayFunc::Fill(weights, (real)1);
		}

		////////////////////////////////////////////////////////////////////////
		//Operations on the global stiffness matrix
		template<int d> void Add_Cell_Stiffness_Matrix(/*rst*/SparseMatrix<real>& K, const MatrixX& K_e, const Array<int>& node_indices)
		{
			const int node_n = (int)node_indices.size(); for (int Ke_i = 0; Ke_i < node_n; Ke_i++) {
				int K_i = node_indices[Ke_i]; for (int Ke_j = 0; Ke_j < node_n; Ke_j++) { int K_j = node_indices[Ke_j]; SparseFunc::Add_Block<d>(K, K_i, K_j, K_e, Ke_i, Ke_j); }
			}
		}
		template void Add_Cell_Stiffness_Matrix<2>(/*rst*/SparseMatrix<real>& K, const MatrixX& K_e, const Array<int>& nodes);
		template void Add_Cell_Stiffness_Matrix<3>(/*rst*/SparseMatrix<real>& K, const MatrixX& K_e, const Array<int>& nodes);

		void Set_Dirichlet_Boundary_Helper(SparseMatrix<real>& K, VectorX& b, const int i, const real psi_D_value)
		{
			for (InnerIterator<real> iter(K, i); iter; ++iter) {
				int j = (int)iter.col();
				if (i == j) { b(j) = psi_D_value * K.coeff(i, j); }
				else {
					real K_ij = K.coeff(i, j); b(j) -= K_ij * psi_D_value;
					K.coeffRef(i, j) = (real)0; K.coeffRef(j, i) = (real)0;
				}
			}
		}
	}
}
