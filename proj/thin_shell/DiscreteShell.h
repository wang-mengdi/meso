//////////////////////////////////////////////////////////////////////////
// Nonlinear Thin Shell FEM
// Copyright (c) (2021-), Fan Feng
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "Hashtable.h"
#include "ElasticParam.h"
#include "BoundaryConditionMesh.h"
#include "Mesh.h"
#include "DiscreteShellParticles.h"
#include "DiscreteShellAuxFunc.h"
#include "Simulator.h"
#include "Common.h"
#include "SparseFunc.h"

using namespace Meso;
enum class AdvanceMode { Explicit, Implicit, Quasistatic };

template<int d> class DiscreteShell : public Meso::Simulator
{
	Typedef_VectorD(d); Typedef_MatrixD(d); 
public:
	AdvanceMode advance_mode=AdvanceMode::Implicit;
	DiscreteShellParticles<d> particles;
	std::shared_ptr< SurfaceMesh<d> > mesh=nullptr;
	HashtableMultiValue<Vector<int,d-1>,int> edge_element_hashtable; //edge (vertices indices: 2 for 3d, 1 for 2d) to face index

	Array<real> hs;					//thickness of the thin shell,defined on vertices
	DiagonalMatrix<real> M;			//mass on vertices
	SparseMatrix<real> A;			//stiffness matrix
	Array<real> dx;					//displacement, delta x
	Array<real> b;					//rhs

	real density;					//density of the thin shell, may need to be moved to materials
	real thickness;					//initial thickness
	Array<ElasticParam> materials;	//different mateirals
	Array<int> material_id;			//refer to the index in materials
	BoundaryConditionMesh<d> bc;	
	Array<MatrixD> Dm_inv;			//[X_j−X_i,X_k−X_i,N]^{-1}
	Array<VectorD> X0;				//positions of rest particles
	Array<real> areas_hat;			//area of faces
	Array<VectorD> N;				//normals of rest triangles
	Array<Vector2i> edges;			//edges in 3d
	Array<real> theta_hats;			//theta at resting position, defined on edges
	Array<real> lambdas;			//coefficient for bending, defined on edges

	//intermediate variables for calculating bending hessian, only for the 3 dimensional case
	Array<ArrayF<Vector3, 3>> ls;
	Array<Vector3> l_norms;
	Array<real> areas;
	Array<Vector3> ps;
	Array<Vector3> hts;
	Array<ArrayF<Vector2,3>> rs; //two indicates the (+-) sign in the formula
	Array<ArrayF<Matrix3, 3>> grad_ns;
	Array<Array2DF<ArrayF<Vector3, 2>,3,3>> grad_rs; //2 grads (+-) of Vector3 for each q and each vertex in a triangle
	
	bool use_exact_hessian = true;		//turn on when the solver is used for optimization, otherwise the we use Gauss-Newton approximation for the forward solve 
	bool use_body_force=true;
	VectorD g=VectorD::Unit(1)*-9.8;
	real damping;						//damping constant

	virtual real CFL_Time(const real cfl);
	virtual void Output(DriverMetaData& metadata);
	virtual void Advance(DriverMetaData& metadata);
	void Initialize(SurfaceMesh<d>& _mesh);
	void Allocate_A();
	void Initialize_Material();

	//Three advancing methods
	void Advance_Explicit(const real dt);
	void Advance_Implicit(const real dt);
	void Advance_Quasistatic();

	void Add_Material(real youngs,real poisson);
	void Set_Fixed(const int node);
	void Set_Displacement(const int node,const VectorD& dis);
	void Set_Force(const int node,const VectorD& force);
	void Add_Force(const int node,const VectorD& force);
	void Clear_Force();
	void Set_Rest_Shape(const Array<VectorD>& _X0);

	//System assembling for the implicit solver, used in Advance_Implicit
	void Clear_A_And_Rhs();
	void Update_Implicit_Force_And_Mass(const real dt=1.0);
	void Update_Implicit_Stretching(const real dt=1.0);
	void Update_Implicit_Bending(const real dt=1.0);
	void Update_Implicit_Boundary_Condition();
	void Solve();

	//Computation of physical quantities
	void Grad_Stretch(const int ele_idx, MatrixD& grad);
	void Grad_Stretch(real area, const MatrixD& stress, const MatrixD&x_hat, const MatrixD& Dm_inv, MatrixD& grad);
	void Grad_Bend(int jct_idx, Eigen::Matrix<real,d,d+1>& grad, const ArrayF<int,d+1>& vtx_idx, const ArrayF<int, 2>& ele_idx);
	void Grad_Bend(const Eigen::Matrix<real,d,d+1>& dtheta_dx, real theta, real theta_hat, real lambda, Eigen::Matrix<real,d,d+1>& dE_dx);
	void Grad_Theta(Eigen::Matrix<real, d, d + 1>& grad_theta, const ArrayF<int, d + 1>& vtx_idx, const ArrayF<int, 2>& ele_idx);
	//Calculate both gradient and hessian
	void Stretch_Force(const real& area, const real& ks, const real& poisson_ratio, const MatrixD& stress, const MatrixD&x_hat, const MatrixD& Dm_inv, MatrixD& grad, MatrixD hess[d][d]);
	void Stretch_Force(const int ele_idx, MatrixD& grad_s, MatrixD hess_s[d][d]);
	void Bend_Force_Approx(int edge_idx, Eigen::Matrix<real,d,d+1>& grad, const ArrayF<int,d+1>& vtx_idx, const ArrayF<int, 2>& ele_idx, MatrixD hess [d+1][d+1]);
	void Bend_Force(int jct_idx, Eigen::Matrix<real, d, d + 1>& grad, const ArrayF<int, d + 1>& vtx_idx, const ArrayF<int, 2>& ele_idx, MatrixD hess[d + 1][d + 1]);
	void Update_Bending_Hess_Variables(); //iterate through elements to update variables used in accurate hessian

	inline real Stretching_Energy(int ele_idx);
	inline real Stretching_Energy(real a, real ks, real nu, const Matrix<real, d>& strain);
	inline real Bending_Energy(real lambda, real theta, real theta_hat);
	inline real Bending_Energy(int jct_idx, const ArrayF<int, d + 1>& vtx_idx, const ArrayF<int, 2>& ele_idx);
	real Total_Stretching_Energy();
	real Total_Bending_Energy();
	real Total_Potential_Energy();
	real Total_Energy();
	inline real Ks(real youngs, real density, real poisson) const { return youngs * density / ((real)1 - poisson * poisson); }
	inline real Kb(real density, real ks) const { return ks * density * density / (real)12; }
	inline real Lambda(real kb, real a_hat, real l_hat = (real)1);
	real Lambda(const ArrayF<int, d + 1>& vtx_idx, const ArrayF<int, 2>& ele_idx);
	void Reset_To_Rest_Position();

	////Helper functions
	const MatrixD C_c () const;
	inline void Strain_To_Stress(real ks, real poisson_ratio, const MatrixD& strain, MatrixD& stress);
	inline void Deformation_To_Strain(const MatrixD& deformation,MatrixD& strain);
	inline real Mass(const int i) const {return M.diagonal()[i*d];}
	inline Array<VectorD>& X(){return particles.xRef();}
	inline Array<VectorD>& V(){return particles.vRef();}
	inline Array<VectorD>& F(){return particles.fRef();}
	inline auto& E(){return mesh->Elements();}
	inline const Array<VectorD>& X() const {return particles.xRef();}
	inline const Array<VectorD>& V() const {return particles.vRef();}
	inline const Array<VectorD>& F() const {return particles.fRef();}
	inline const auto& E() const {return mesh->Elements();}
	inline int Vtx_Num() const {return particles.Size();}
	inline int Ele_Num() const {return (int)mesh->Elements().size();}
	inline int Jct_Num() const { if constexpr (d == 2) { return Vtx_Num(); } else { return edges.size(); } }	//junction size
	bool Junction_Info(int junction_idx, ArrayF<int, d + 1>& vtx_idx, ArrayF<int,2>& ele_idx);  //junction index is the edge index in 3d and the vertex indx in 2d, where the dihedral angle is defined
	inline real Ele_H(int ele_idx) const { real avg_h = 0; for (int i = 0; i < d; i++) { avg_h += hs[E()[ele_idx][i]]; } return avg_h/(real)d; } //avg thickness of an element
	inline real Jct_H(int jct_idx) const {if constexpr (d == 3) {return (real)0.5 * (hs[edges[jct_idx][0]] + hs[edges[jct_idx][1]]);}else { return hs[jct_idx]; }	} //thickness of jct

	//Numerical Verification
	MatrixD Numerical_Grad_Stretch(int ele_idx,const MatrixD& grad_s);
	Array2DF<Matrix<real, d>, d, d> Numerical_Hess_Stretch(int ele_idx, const MatrixD& grad_s, const MatrixD hess_s[d][d]);
	Eigen::Matrix<real,d,d+1> Numerical_Grad_Bend(int jct_idx, ArrayF<int, d + 1>& vtx_idx, ArrayF<int, 2>& ele_idx, const Eigen::Matrix<real, d, d + 1>& grad_b);
	Array2DF<MatrixD, d + 1, d + 1> Numerical_Hess_Bend(int jct_idx, const ArrayF<int, d + 1>& vtx_idx,const ArrayF<int, 2>& ele_idx, const Eigen::Matrix<real,d,d+1>& grad_b, const MatrixD hess_b[d + 1][d + 1]);
	Array2DF<MatrixD, d + 1, d + 1> Numerical_Hess_Theta(const ArrayF<int, d + 1>& vtx_idx,const ArrayF<int, 2>& ele_idx, const Eigen::Matrix<real, d, d + 1>& grad_theta);

protected:
	void Add_Block(Array<real>& b, const int i, const VectorD& bi);
};
