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
			dN0/ds dN1/ds dN2/ds dN3/ds
			dN0/dt dN1/dt dN2/dt dN3/dt
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
			//N0=1/8(1-s)(1-t)(1-q)
			//N1=1/8(1+s)(1-t)(1-q)
			//N2=1/8(1+s)(1+t)(1-q)
			//N3=1/8(1-s)(1+t)(1-q)
			//N4=1/8(1-s)(1-t)(1+q)
			//N5=1/8(1+s)(1-t)(1+q)
			//N6=1/8(1+s)(1+t)(1+q)
			//N7=1/8(1-s)(1+t)(1+q)

			/*
			dN0/ds dN1/ds dN2/ds dN3/ds dN4/ds dN5/ds dN6/ds dN7/ds
			dN0/dt dN1/dt dN2/dt dN3/dt dN4/dt dN5/dt dN6/dt dN7/dt
			dN0/dq dN1/dq dN2/dq dN3/dq dN4/dq dN5/dq dN6/dq dN7/dq
			*/

			/*			far		t  q
			   N7 --  N6		| /
			  /|	 /|			|/
			 / |    / |			-->s
			N3-|- N2  |
			|  N4 |- N5
			| /	  |  /
			|/	  | /
			N0 -- N1
			close
			*/

			const real s = natural_coord[0]; const real t = natural_coord[1]; const real q = natural_coord[2];
			dNde << -(1 - t) * (1 - q), (1 - t) * (1 - q), (1 + t) * (1 - q), -(1 + t) * (1 - q), -(1 - t)* (1 + q), (1 - t)* (1 + q), (1 + t)* (1 + q), -(1 + t) * (1 + q),
					-(1 - s) * (1 - q), -(1 + s) * (1 - q), (1 + s)* (1 - q), (1 - s)* (1 - q), -(1 - s) * (1 + q), -(1 + s) * (1 + q), (1 + s)* (1 + q), (1 - s)* (1 + q),
					-(1 - s) * (1 - t), -(1 + s)* (1 - t), -(1 + s)* (1 + t), -(1 - s) * (1 + t), (1 - s) * (1 - t), (1 + s)* (1 - t), (1 + s)* (1 + t), (1 - s)* (1 + t);
			dNde *= (real).125;
		}

		template<int d> void Cell_dNdX(const Vector<real, d>& natural_coord, const real dx, MatrixX& dNdX)
		{
			dNde<d>(natural_coord, dNdX); dNdX *= ((real)2 / dx); //dn/dx=dn/de*de/dx, de/dx=2/dx
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
			/*
				dN0/dx	0	0	dN1/dx	0	0	dN2/dx	0	0	dN3/dx	0	0	dN4/dx	0	0	dN5/dx	0	0	dN6/dx	0	0	dN7/dx	0	0
				0	dN0/dy	0	0	dN1/dy	0	0	dN2/dy	0	0	dN3/dy	0	0	dN4/dy	0	0	dN5/dy	0	0	dN6/dy	0	0	dN7/dy	0
				0	0	dN0/dz	0	0	dN1/dz	0	0	dN2/dz	0	0	dN3/dz	0	0	dN4/dz	0	0	dN5/dz	0	0	dN6/dz	0	0	dN7/dz
				dN0/dy	dN0/dx	0	dN1/dy	dN1/dx	0	N2/dy	dN2/dx	0	dN3/dy	dN3/dx	0	dN4/dy	dN4/dx	0	dN5/dy	dN5/dx	0	N6/dy	dN6/dx	0	dN7/dy	dN7/dx	0
				0	dN0/dz	dN0/dy	0	dN1/dz	dN1/dy	0	N2/dz	dN2/dy	0	dN3/dz	dN3/dy	0	dN4/dz	dN4/dy	0	dN5/dz	dN5/dy	0	N6/dz	dN6/dy	0	dN7/dz	dN7/dy
				dN0/dz	0	dN0/dx	dN1/dz	0	dN1/dx	N2/dz	0	dN2/dx	dN3/dz	0	dN3/dx	dN4/dz	0	dN4/dx	dN5/dz	0	dN5/dx	N6/dz	0	dN6/dx	dN7/dz	0	dN7/dx
			*/
			int vtx_n = (int)pow(2, 3);
			B.resize(6, 24); B.fill(0);
			MatrixX dNdX(3, vtx_n); Cell_dNdX<3>(natural_coord, dx, dNdX);

			const int x = 0; const int y = 1; const int z = 2;
			for (int i = 0; i < 8; i++) {
				B(0, i * 3) = dNdX(x, i);
				B(1, i * 3 + 1) = dNdX(y, i);
				B(2, i * 3 + 2) = dNdX(z, i);
				B(3, i * 3) = dNdX(y, i); B(3, i * 3+1) = dNdX(x, i);
				B(4, i * 3 + 1) = dNdX(z, i); B(4, i * 3 + 2) = dNdX(y, i);
				B(5, i * 3 ) = dNdX(z, i);	B(5, i * 3+2) = dNdX(x, i);
			}
		}

		template<int d> void Cell_Stiffness_Matrix(const real youngs, const real poisson, const real dx, MatrixX& K_e)
		{
			int n = d * (int)pow(2, d); K_e.resize(n, n); K_e.fill((real)0);
			Matrix<real,1<<d,d> points; Vector<real,1<<d> weights; 
			Initialize_Gaussian_Integration_Points<d>(points, weights);
			MatrixX E; Strain_Stress_Matrix_Linear<d>(youngs, poisson, E);
			real J_det = pow(dx / (real)2, (real)d);
			for (auto i = 0; i < points.rows(); i++) {
				MatrixX B; Cell_Strain_Displacement_Matrix<d>(points.row(i), dx, B);
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
			K_e.resize(24, 24); K_e.fill((real)0);
			VectorX k(14);
			k << -((real)6*poisson -(real) 4) /9, (real)1 / (real)12, -(real)1 /(real) 9, -((real)4 * poisson - (real)1)/ (real)12, ((real)4 * poisson -(real)1) / 12, (real)1 / (real)18, (real)1 / (real)24,
				-(real)1 / (real)12, ((real)6* poisson -(real)5) / (real)36, -((real)4* poisson - (real)1) / (real)24, -(real)1 / (real)24, ((real)4* poisson-(real)1) / (real)24, ((real)3* poisson -(real)1) / (real)18, ((real)3* poisson -(real) 2) / (real)18;
			MatrixX k0(6, 6); MatrixX k1(6, 6); MatrixX k2(6, 6); MatrixX k3(6, 6); MatrixX k4(6, 6); MatrixX k5(6, 6);
			/*k1 k2 k2 k3 k5 k5
			  k2 k2 k1 k4 k7 k6
			  k3 k4 k4 k1 k8 k8
			  k5 k6 k7 k8 k1 k2
			  k5 k7 k6 k8 k2 k1*/
			k0.row(0) << k[0], k[1] ,k[1] ,k[2] ,k[4] , k[4];
			k0.row(1) << k[1], k[0], k[1], k[3], k[5], k[6];
			k0.row(2) << k[1], k[1], k[0], k[3], k[6], k[5];
			k0.row(3) << k[2], k[3], k[3], k[0], k[7], k[7];
			k0.row(4) << k[4], k[5], k[6], k[7], k[0], k[1];
			k0.row(5) << k[4], k[6], k[5], k[7], k[1], k[0];
			/*k9 k8 k12 k6 k4 k7
			  k8 k9 k12 k5 k3 k5
			  k10 k10 k13 k7 k4 k6
			  k6 k5 k11 k9 k2 k10
			  k4 k3 k5 k2 k9 k12
			  k11 k4 k6 k12 k10 k13*/
			k1.row(0) << k[8], k[7], k[11], k[5], k[3], k[6];
			k1.row(1) << k[7], k[8], k[11], k[4], k[2], k[4];
			k1.row(2) << k[9], k[9], k[12], k[6], k[3], k[5];
			k1.row(3) << k[5], k[4], k[10], k[8], k[1], k[9];
			k1.row(4) << k[3], k[2], k[4], k[1], k[8], k[11];
			k1.row(5) << k[10], k[3], k[5], k[11], k[9], k[12];
			/*k6 k7 k4 k9 k12 k8
			  k7 k6 k4 k10 k13 k10
			  k5 k5 k3 k8 k12 k9
			  k9 k10 k2 k6 k11 k5
			  k12 k13 k10 k11 k6 k4
			  k2 k12 k9 k4 k5 k3*/
			k2.row(0) << k[5], k[6], k[3], k[8], k[11], k[7];
			k2.row(1) << k[6], k[5], k[3], k[9], k[12], k[9];
			k2.row(2) << k[4], k[4], k[2], k[7], k[11], k[8];
			k2.row(3) << k[8], k[9], k[1], k[5], k[10], k[4];
			k2.row(4) << k[11], k[12], k[9], k[10], k[5], k[3];
			k2.row(5) << k[1], k[11], k[8], k[3], k[4], k[2];
			/*k14 k11 k11 k13 k10 k10
			  k11 k14 k11 k12 k9 k8
			  k11 k11 k14 k12 k8 k9
			  k13 k12 k12 k14 k7 k7
			  k10 k9 k8 k7 k14 k11
			  k10 k8 k9 k7 k11 k14*/
			k3.row(0) << k[13], k[10], k[10], k[12], k[9], k[9];
			k3.row(1) << k[10], k[13], k[10], k[11], k[8], k[7];
			k3.row(2) << k[10], k[10], k[13], k[11], k[7], k[8];
			k3.row(3) << k[12], k[11], k[11], k[13], k[6], k[6];
			k3.row(4) << k[9], k[8], k[7], k[6], k[13], k[10];
			k3.row(5) << k[9], k[7], k[8], k[6], k[10], k[13];
			/*k1 k2 k8 k3 k5 k4
			  k2 k1 k8 k4 k6 k11
			  k8 k8 k1 k5 k11 k6
			  k3 k4 k5 k1 k8 k2
			  k5 k6 k11 k8 k1 k8
			  k4 k11 k6 k2 k8 k1*/
			k4.row(0) << k[0], k[1], k[7], k[2], k[4], k[3];
			k4.row(1) << k[1], k[0], k[7], k[3], k[5], k[10];
			k4.row(2) << k[7], k[7], k[0], k[4], k[10], k[5];
			k4.row(3) << k[2], k[3], k[4], k[0], k[7], k[1];
			k4.row(4) << k[4], k[5], k[10], k[7], k[0], k[7];
			k4.row(5) << k[3], k[10], k[5], k[1], k[7], k[0];
			/*k14 k11 k7 k13 k10 k12
			  k11 k14 k7 k12 k9 k2
			  k7 k7 k14 k10 k2 k9
			  k13 k12 k10 k14 k7 k11
			  k10 k9 k2 k7 k14 k7
			  k12 k2 k9 k11 k7 k14*/
			k5.row(0) << k[13], k[10], k[6], k[12], k[9], k[11];
			k5.row(1) << k[10], k[13], k[6], k[11], k[8], k[1];
			k5.row(2) << k[6], k[6], k[13], k[9], k[1], k[8];
			k5.row(3) << k[12], k[11], k[9], k[13], k[6], k[10];
			k5.row(4) << k[9], k[8], k[1], k[6], k[13], k[6];
			k5.row(5) << k[11], k[1], k[8], k[10], k[6], k[13];
			
			/*k1 k2 k3 k4
			kT2 k5 k6 kT4
			kT3 k6 kT5 kT2
			k4 k3 k2 kT1*/
			K_e.block<6, 6>(0, 0) = k0; K_e.block<6, 6>(0, 6) = k1;	K_e.block<6, 6>(0, 12) = k2; K_e.block<6, 6>(0, 18) = k3;
			K_e.block<6, 6>(6, 0) = k1.transpose(); K_e.block<6, 6>(6, 6) = k4;	K_e.block<6, 6>(6, 12) = k5; K_e.block<6, 6>(6, 18) = k2.transpose();
			K_e.block<6, 6>(12, 0) = k2.transpose(); K_e.block<6, 6>(12, 6) = k5;	K_e.block<6, 6>(12, 12) = k4.transpose(); K_e.block<6, 6>(12, 18) = k1.transpose();
			K_e.block<6, 6>(18, 0) = k3; K_e.block<6, 6>(18, 6) = k2;	K_e.block<6, 6>(18, 12) = k1; K_e.block<6, 6>(18, 18) = k0.transpose();

			K_e *= youngs;
			K_e /= (((real)1 +poisson)*((real)1 - (real)2*poisson));
			K_e /= (real)2; //account for the change of coordinate
		}

		////////////////////////////////////////////////////////////////////////
		//Gaussian integration
		template<> void Initialize_Gaussian_Integration_Points<2>(Matrix<real, 4, 2>& points, Vector<real,4>& weights) //2^2=4 sample points
		{
			real c = (real)1 / sqrt((real)3); points.row(0) = Vector2(-c, -c); points.row(1) = Vector2(c, -c); points.row(2) = Vector2(-c, c); points.row(3) = Vector2(c, c); weights.fill(1);
		}

		template<> void Initialize_Gaussian_Integration_Points<3>(Matrix<real, 8, 3>& points, Vector<real, 8>& weights) //2^3=8 sample points
		{
			real c = (real)1 / sqrt((real)3); points.row(0) = Vector3(-c, -c, -c); points.row(1) = Vector3(c, -c, -c); points.row(2) = Vector3(c, c, -c); points.row(3) = Vector3(-c, c, -c);
			points.row(4) = Vector3(-c, -c, c); points.row(5) = Vector3(c, -c, c); points.row(6) = Vector3(c, c, c); points.row(7) = Vector3(-c, c, c); weights.fill(1);
		}

		////////////////////////////////////////////////////////////////////////
		//Operations on the global stiffness matrix
		template<int d> void Add_Cell_Stiffness_Matrix(/*rst*/SparseMatrix<real>& K, const MatrixX& K_e, const Array<int>& node_indices)
		{
			const int node_n = (int)node_indices.size(); 
			for (int Ke_i = 0; Ke_i < node_n; Ke_i++) {
				int K_i = node_indices[Ke_i]; 
				for (int Ke_j = 0; Ke_j < node_n; Ke_j++) { 
					int K_j = node_indices[Ke_j];
					SparseFunc::Add_Block<real,d,d>(K, K_i, K_j, K_e, Ke_i, Ke_j); 
				}
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
