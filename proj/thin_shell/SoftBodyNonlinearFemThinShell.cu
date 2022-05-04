//////////////////////////////////////////////////////////////////////////
// Nonlinear Thin Shell FEM
// Copyright (c) (2021-), Fan Feng
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#include "SoftBodyNonlinearFemThinShell.h"
#include "MeshFunc.h"
#include "AuxFunc.h"
#include "NonlinearFemFunc.h"
#include "Timer.h"
#include "Hashtable.h"
#include "SimplicialPrimitives.h"
#include "IOHelper.h"
#include "Common.h"
#include <Eigen/IterativeLinearSolvers>
#include <iostream>
#include <fstream>

using namespace ThinShellAuxFunc; 
using namespace Meso;

template<class T_ARRAY> int Element_Edges(const Vector2i& v,T_ARRAY& edges);
template<class T_ARRAY> int Element_Edges(const Vector3i& v,T_ARRAY& edges);
void Grad_Q(const ArrayF<Vector3, 3>& vtx, const int i, const int j, const Vector3& ps, const real& qs_i, const ArrayF<Vector3, 3>& ls, const real& a, Vector3& grad_q);
void Grad_R(const ArrayF<Vector3, 3>& vtx, const int i, const int j, const Vector3& ps, const ArrayF<Vector3, 3>& ls, const real& a, ArrayF<Vector3, 2>& grad_r);
void Grad_N(const ArrayF<Vector3, 3>& vtx, const ArrayF<Vector3, 3>& ls, const real& a, ArrayF<Matrix3, 3>& grad_n);

template<int d> real SoftBodyNonlinearFemThinShell<d>::CFL_Time(const real cfl) {
	return 0.01;
}

template<int d> void SoftBodyNonlinearFemThinShell<d>::Output(const bf::path base_path, const int frame) {
	std::string vtu_name = fmt::format("vtu{:04d}.vtu", frame);
	bf::path vtu_path = base_path / bf::path(vtu_name);
	VTKFunc::Output_VTU<d, VectorD>(mesh, V(), vtu_path.string());
}

template<int d> void SoftBodyNonlinearFemThinShell<d>::Advance(const int current_frame, const real current_time, const real dt) {
	if (use_explicit) { Advance_Explicit(dt); }
	else { Advance_Implicit(dt); }
	return;
}

template<int d> void SoftBodyNonlinearFemThinShell<d>::Initialize(SurfaceMesh<d>& _mesh)
{
	mesh=std::make_shared<SurfaceMesh<d>>(particles.XPtr());
	particles.Resize((int)_mesh.Vertices().size());
	*mesh=_mesh;
	MeshFunc::Get_Edges<d,d>(_mesh,edges); //no repetition in edges
	
	//initialize edge hashtable
	for(int i=0;i<Ele_Num();i++){
		VectorDi& vtx_indices=E()[i];
		ArrayF<Vector<int,d-1>,d> edge_set; 
		Element_Edges(vtx_indices,edge_set);
		for(int k=0;k<d;k++){
			Add(edge_element_hashtable,Unique_Ordered(edge_set[k]),i);
		}
	}

	int vtx_n=Vtx_Num();int dof_n=vtx_n*d;int ele_n=Ele_Num();
	material_id.resize(ele_n,0);

	////initialize X0, Dm_inv, N and M
	Set_Rest_Shape(particles.XRef());
	M.resize(dof_n);for(int i=0;i<M.diagonal().size();i++)M.diagonal()[i]=(real)0;
	Dm_inv.resize(ele_n,MatrixD::Zero());
	N.resize(ele_n, VectorD::Zero());
	areas_hat.resize(ele_n);

	ArrayFunc::Fill(F(),VectorD::Zero());
	ArrayFunc::Fill(V(),VectorD::Zero());

	////initialize implicit variables
	if (!use_explicit) {
		A.resize(dof_n,dof_n);
		Allocate_A();
		dv.resize(dof_n);dv.fill((real)0);
		b.resize(dof_n);b.fill((real)0);

		if constexpr (d == 3) {
			if (use_exact_hessian) {
				ls.resize(ele_n);//fill?
				l_norms.resize(ele_n);
				areas.resize(ele_n);
				ps.resize(ele_n, Vector3::Zero());
				hts.resize(ele_n, Vector3::Zero());
				rs.resize(ele_n);//fill?
				grad_ns.resize(ele_n);
				grad_rs.resize(ele_n);
			}
		}
	}
}

template<int d> void SoftBodyNonlinearFemThinShell<d>::Allocate_A() {
	std::vector<Triplet<real>> triplets;

	//vertex with itself
	for (int i = 0; i < Vtx_Num(); i++) {
		int r = i; int c = i;
		for (int rr = r * d; rr < (r + 1) * d; rr++)for (int cc = c * d; cc < (c + 1) * d; cc++) {
			triplets.push_back(Triplet<real>(rr, cc, (real)0));
		}
	}

	//vertex with neighboring vertices
	for (int i = 0; i < edges.size(); i++) {
		const Vector2i& e = edges[i]; int r = e[0]; int c = e[1];
		for (int rr = r * d; rr < (r + 1) * d; rr++)for (int cc = c * d; cc < (c + 1) * d; cc++) {
			triplets.push_back(Triplet<real>(rr, cc, (real)0));
		}
		r = e[1]; c = e[0];
		for (int rr = r * d; rr < (r + 1) * d; rr++)for (int cc = c * d; cc < (c + 1) * d; cc++) {
			triplets.push_back(Triplet<real>(rr, cc, (real)0));
		}

		if constexpr (d == 3) {
			Array<int> incident_elements;
			Value_Array(edge_element_hashtable, edges[i], incident_elements);
			if (incident_elements.size() == 2) {
				int face_idx_0 = (incident_elements[0]), face_idx_1 = (incident_elements[1]);
				
				//Find the two vertices at opposite sides
				r = Third_Vertex(e[0], e[1], E()[face_idx_0]);
				c = Third_Vertex(e[0], e[1], E()[face_idx_1]);
				Assert(r != -1 && c != -1, "index {} {} out of range for finding opposite vertex");
				for (int rr = r * d; rr < (r + 1) * d; rr++)for (int cc = c * d; cc < (c + 1) * d; cc++) {
					triplets.push_back(Triplet<real>(rr, cc, (real)0));
					triplets.push_back(Triplet<real>(cc, rr, (real)0));
				}
			}
		}
	}

	A.setFromTriplets(triplets.begin(), triplets.end());
}

template<int d> void SoftBodyNonlinearFemThinShell<d>::Initialize_Material() {
	//initialize thickness
	hs.resize(Vtx_Num(), thickness);
	for (int i = 0; i < Ele_Num(); i++) {
		ArrayF<VectorD, d> vtx;
		for (int j = 0; j < d; j++) { vtx[j] = X0[E()[i][j]]; }
		NonlinearFemFunc<d>::D_Inv_And_Area_And_Normal(vtx, Dm_inv[i], areas_hat[i], N[i]);

		for (int j = 0; j < d; j++) {
			int v_idx = E()[i][j];
			for (int k = 0; k < d; k++) {
				M.diagonal()[v_idx * d + k] += hs[v_idx] * areas_hat[i] / (real)(d);	//Mass on vertices
			}
		}
	}

	//initialize lambda, theta
	if constexpr (d == 3) {
		lambdas.resize(edges.size(),(real)1);
		theta_hats.resize(edges.size(),(real)0);

		for(int edge_idx=0; edge_idx <edges.size(); edge_idx++){
			ArrayF<int,d+1> vtx_idx;
			ArrayF<int, 2> ele_idx;
			if (Junction_Info(edge_idx, vtx_idx, ele_idx)) { //shared edge
				real l_hat = (X()[vtx_idx[0]] - X()[vtx_idx[1]]).norm();
				real a_hat = areas_hat[ele_idx[0]]+ areas_hat[ele_idx[1]];
				ElasticParam& material0=materials[material_id[ele_idx[0]]];
				ElasticParam& material1=materials[material_id[ele_idx[1]]];
				real avg_h0, avg_h1, avg_h;
				avg_h0 = avg_h1 = avg_h = 0;
				for (int i = 0; i < 4; i++) { avg_h += hs[vtx_idx[i]]; } avg_h /= (real)4;
				for (int i = 0; i < 3; i++) { avg_h0 += hs[E()[ele_idx[0]][i]]; } avg_h0 /= (real)3;
				for (int i = 0; i < 3; i++) { avg_h1 += hs[E()[ele_idx[1]][i]]; } avg_h1 /= (real)3;
				real ks0=Ks(material0.youngs_modulus, avg_h,material0.poisson_ratio);	//this thickness should use the thickness on face instead?
				real ks1=Ks(material1.youngs_modulus, avg_h,material1.poisson_ratio);	//although it really simplifies the calculation
				lambdas[edge_idx]=Lambda(Kb(avg_h,(real)0.5*(ks0+ks1)),a_hat,l_hat);
				VectorD n0=Triangle<3>::Normal(X()[E()[ele_idx[0]][0]],X()[E()[ele_idx[0]][1]],X()[E()[ele_idx[0]][2]]);
				VectorD n1=Triangle<3>::Normal(X()[E()[ele_idx[1]][0]],X()[E()[ele_idx[1]][1]],X()[E()[ele_idx[1]][2]]);
				theta_hats[edge_idx]=Dihedral_Angle(n0,n1,X()[vtx_idx[0]],X()[vtx_idx[1]]);
			}
		}
	}
}

template<int d> void SoftBodyNonlinearFemThinShell<d>::Add_Material(real youngs,real poisson)
{materials.push_back(ElasticParam(youngs,poisson));}

template<int d> void SoftBodyNonlinearFemThinShell<d>::Set_Fixed(const int node)
{bc.psi_D_values[node]=VectorD::Zero();}

template<int d> void SoftBodyNonlinearFemThinShell<d>::Set_Displacement(const int node,const VectorD& dis)
{bc.psi_D_values[node]=dis;}

template<int d> void SoftBodyNonlinearFemThinShell<d>::Set_Force(const int node,const VectorD& force)
{bc.forces[node]=force;}

template<int d> void SoftBodyNonlinearFemThinShell<d>::Add_Force(const int node,const VectorD& force)
{bc.forces[node]+=force;}

template<int d> void SoftBodyNonlinearFemThinShell<d>::Clear_Force()
{bc.forces.clear();}

template<int d> void SoftBodyNonlinearFemThinShell<d>::Set_Rest_Shape(const Array<VectorD>& _X0)
{X0=_X0;}

template<int d> void SoftBodyNonlinearFemThinShell<d>::Advance_Explicit(const real dt)
{
	Timer timer;

	ArrayFunc::Fill(F(),VectorD::Zero());
	const int vtx_num=Vtx_Num();

	////body force, damping force, boundary
	if(use_body_force){for(int i=0;i<vtx_num;i++)F()[i]+=Mass(i)*g; }

	//damping
	for(int i=0;i<vtx_num;i++)F()[i]-=Mass(i)*damping*V()[i];

	//stretching forces
	const int ele_num = Ele_Num();
	for(int ele_idx=0; ele_idx < ele_num; ele_idx++){
		MatrixD grad_s;
		Grad_Stretch(ele_idx,grad_s);
		for(int j=0;j<d;j++){F()[E()[ele_idx][j]]-=grad_s.col(j);}

		/*real stretching_energy=Stretching_Energy<d>(ele_idx);
		MatrixD grad_s_n=Numerical_Grad_Stretch(ele_idx,stretching_energy,grad_s);*/
	}

	////bending forces, only support for three dimension for now
	if constexpr (d == 3) {
		for(int i=0;i<edges.size();i++){
			Eigen::Matrix<real,d,d+1> grad_b;
			ArrayF<int,d+1> vtx_idx;
			ArrayF<int, 2> ele_idx;
			if (Junction_Info(i, vtx_idx, ele_idx)) {
				Grad_Bend(i, grad_b, vtx_idx, ele_idx);
				for(int j=0;j<d+1;j++){
					F()[vtx_idx[j]]-=grad_b.col(j);
				}

				/*Eigen::Matrix<real,d,d+1> grad_b_n=Numerical_Grad_Bend(i,theta_hats[i],lambdas[i]);
				std::cout<<"grad_b["<<i<<"]: \n"<<grad_b<<"\n"<<"grad_b_n["<<i<<"]: \n"<<grad_b_n<<std::endl;*/
			}
		}
	}

	//enforce boundary condition
	for(auto& iter:bc.forces){F()[iter.first]+=iter.second;}

	////time integration
	for(int i=0;i<vtx_num;i++){
		V()[i]+=F()[i]/Mass(i)*dt;
		if(bc.Is_Psi_D(i))V()[i]=VectorD::Zero();
		X()[i]+=V()[i]*dt;}
	Info("Explicit time integration: {} ms" , timer.Total_Time(PhysicalUnits::ms));
}

// A = dt^2 J + dt damp J
// b = dt f + dt^2 J v + dt damp J v
template<int d> void SoftBodyNonlinearFemThinShell<d>::Advance_Implicit(const real dt)
{
	Timer timer;

	Clear_A_And_Rhs();
	std::cout << "# number of non-zeros of A:  " << A.nonZeros() << std::endl;
	const int vtx_num=Vtx_Num();
	const int ele_num=Ele_Num();

	Update_Implicit_Force_And_Mass(dt);
	Update_Implicit_Stretching(dt);
	Update_Implicit_Bending(dt);
	Update_Implicit_Boundary_Condition(dt);

	Eigen::ConjugateGradient<SparseMatrix<real>, Eigen::Lower | Eigen::Upper, Eigen::DiagonalPreconditioner<real>> cg;
	cg.setTolerance((real)1e-5);
	dv = cg.compute(A).solve(b);
	Info("Linear system solve: {} ms", timer.Lap_Time(PhysicalUnits::ms));

	std::cout << "#	CG iterations:     " << cg.iterations() << std::endl;
	std::cout << "#	CG estimated error: " << cg.error() << std::endl;

	//Eigen::LLT<Eigen::MatrixXd> lltOfA(A); // compute the Cholesky decomposition of A
	//if (lltOfA.info() == Eigen::NumericalIssue){throw std::runtime_error("Negative matrix!");}
	//if (!A.isApprox(A.transpose())) {throw std::runtime_error("Non symmetric matrix!");}

#pragma omp parallel for
	for (int i = 0; i < vtx_num; i++) {
		for (int j = 0; j < d; j++) { V()[i][j] += dv[i * d + j]; }
		X()[i] += V()[i] * dt;
	}

	Info("update nodes: {} ms", timer.Lap_Time(PhysicalUnits::ms));
}

// A = hess
// b = -grad
// Adx=b
template<int d> void SoftBodyNonlinearFemThinShell<d>::Advance_Quasi_Static()
{	//dv is dx in this case
	int iter = 0;
	real err = 1;
	const int vtx_num = Vtx_Num();
	const int ele_num = Ele_Num();
	const real alpha = 0.1;
	real energy = 0;
	const int max_iter = 1000;
	while (err > 1e-4) {
		if (iter == max_iter) {
			Info("max iteration {} reached!",max_iter);
			break;
		}
		//energy = (real)0;
		//Info("");
		//Info("Start the {}th Newton iteration: ", iter);
		//Timer timer;
		//timer.Reset();

		SparseFunc::Set_Value(A, (real)0);
		b.setZero(); //dense vector b can be directly set to zero
		//std::cout << "# number of non-zeros of A:  " << A.nonZeros() << std::endl;

		//add external forces
		for (auto& force : bc.forces) {
			Add_Block(b, force.first, force.second);
			//energy += -force.second.dot(X()[force.first]); //potential energy by the external force
		}
		//timer.Elapse_And_Output_And_Reset("Add external force for b");

		//Stretching
		//Timer timer2;
		//timer2.Begin_Loop();
		for (int ele_idx = 0; ele_idx < ele_num; ele_idx++) {
			MatrixD grad_s; MatrixD hess_s[d][d];
			Stretch_Force(ele_idx, grad_s, hess_s);

			/*MatrixD grad_s_n = Numerical_Grad_Stretch(ele_idx,grad_s);
			Array2DF<Matrix<real, d>, d, d> hess_s_n = Numerical_Hess_Stretch(ele_idx, grad_s,hess_s);*/

			//timer2.Record("calculate streching force");

			//iterate through verteices in the element
			for (int j = 0; j < d; j++) {
				Add_Block(b, E()[ele_idx][j], -grad_s.col(j));
				for (int k = 0; k < d; k++) {
					Add_Block_Helper(A, E()[ele_idx][j], E()[ele_idx][k], hess_s[j][k]);
				}
			}
			//timer2.Record("assemble strething matrix");

			//energy += Stretching_Energy(areas_hat[ele_idx], ks, material.poisson_ratio,strain);
		}

		//timer2.End_Loop_And_Output(std::cout);
		//timer.Elapse_And_Output_And_Reset("Assemble linear system for stretching");

		//Bending
		//timer2.Begin_Loop();
		if constexpr (d == 3) {
			Update_Bending_Hess_Variables();
			for (int jct_idx = 0; jct_idx < edges.size(); jct_idx++) {
				Eigen::Matrix<real, d, d + 1> grad_b;
				MatrixD hess_b[d + 1][d + 1];
				ArrayF<int, d + 1> vtx_idx;
				ArrayF<int, 2> ele_idx;
				if (Junction_Info(jct_idx, vtx_idx, ele_idx)) {
					Bend_Force(jct_idx, grad_b, vtx_idx, ele_idx, hess_b);
					//timer2.Record("calculate bending force");
					for (int j = 0; j < d + 1; j++) {
						Add_Block(b, vtx_idx[j], -grad_b.col(j));
						for (int k = 0; k < d + 1; k++) {
							Add_Block_Helper(A, vtx_idx[j], vtx_idx[k], hess_b[j][k]);
						}
					}

					/*Eigen::Matrix<real,d,d+1> grad_b_n=Numerical_Grad_Bend(jct_idx,vtx_idx,ele_idx,grad_b);
					Array2DF<Matrix<real, d>, d + 1, d + 1> hess_b_n = Numerical_Hess_Bend(jct_idx, vtx_idx, ele_idx, grad_b,hess_b);*/

					//timer2.Record("assemble bending matrix");

					//extra calculation for bending energy, could be simplified with the bend force calculation?
					/*Vector3 n0 = Triangle<3>::Normal(X()[E()[ele_idx[0]][0]], X()[E()[ele_idx[0]][1]], X()[E()[ele_idx[0]][2]]);
					Vector3 n1 = Triangle<3>::Normal(X()[E()[ele_idx[1]][0]], X()[E()[ele_idx[1]][1]], X()[E()[ele_idx[1]][2]]);
					real theta = Dihedral_Angle(n0, n1, X()[vtx_idx[0]], X()[vtx_idx[1]]);

					energy += Bending_Energy(lambdas[jct_idx], theta, theta_hats[jct_idx]);*/
				}
			}
		}
		//timer2.End_Loop_And_Output(std::cout);
		//timer.Elapse_And_Output_And_Reset("Assemble linear system for bending");
		
		for (auto& bc_d : bc.psi_D_values) {
			int node = bc_d.first; VectorD dis = bc_d.second;
			for (int axis = 0; axis < d; axis++) {
				int idx = node * d + axis;
				if(iter==0){ NonlinearFemFunc<d>::Set_Dirichlet_Boundary_Helper(A, b, idx, dis[axis]); }
				else{ NonlinearFemFunc<d>::Set_Dirichlet_Boundary_Helper(A, b, idx, (real)0); }
			}
		}
		
		Eigen::ConjugateGradient<SparseMatrix<real>, Eigen::Lower | Eigen::Upper> cg;
		cg.setTolerance((real)1e-6);
		cg.compute(A);
		dv = cg.solve(b);

		///*std::cout << "A:" << std::endl;
		//std::cout << A << std::endl;
		//std::cout << "b:" << std::endl;
		//std::cout << b.transpose() << std::endl;
		//std::cout << "dv:" << std::endl;
		//std::cout << dv.transpose() << std::endl;*/
		//
		////timer.Elapse_And_Output_And_Reset("linear system solve");

		////std::cout << "#	CG iterations:     " << cg.iterations() << std::endl;
		////std::cout << "#	CG estimated error: " << cg.error() << std::endl;

		//#pragma omp parallel for
		//for (int i = 0; i < vtx_num; i++) {
		//	for (int j = 0; j < d; j++) {
		//		X()[i][j] += damping*dv[i * d + j]; //damping becomes the step size here
		//	}
		//}

		////timer.Elapse_And_Output_And_Reset("update nodes");

		//err = dv.norm() / dv.size();
		////Info("b norm is:{}", b.norm());
		////Info("relative error is:{}", err);
		//iter++;
		////energies_n.push_back(energy);
	}
	Info("Quasi_static solve finished with {} Newton iterations: ", iter);
}

template<int d> void SoftBodyNonlinearFemThinShell<d>::Clear_A_And_Rhs(){
	SparseFunc::Set_Value(A, (real)0);
	b.setZero();
}

template<int d> void SoftBodyNonlinearFemThinShell<d>::Update_Implicit_Force_And_Mass(const real dt){
	Timer timer;
	const int vtx_num=Vtx_Num();
	const int ele_num=Ele_Num();

	//add external forces
	if(use_body_force){
		for (int i = 0; i < vtx_num; i++) {
			Add_Block(b, i, dt * Mass(i) * g);
		}
	}
	Info("Add body force for b", timer.Lap_Time(PhysicalUnits::ms));

	for(auto& iter:bc.forces){
		Add_Block(b,iter.first, dt*iter.second);
	}
	Info("Add external force for b", timer.Lap_Time(PhysicalUnits::ms));

	for(int i=0;i<vtx_num;i++){Add_Block_Helper(A,i,i,Mass(i)*MatrixD::Identity());}
	Info("Add mass for A", timer.Lap_Time(PhysicalUnits::ms));
}

template<int d> void SoftBodyNonlinearFemThinShell<d>::Update_Implicit_Stretching(const real dt){
	Timer timer;
	const int vtx_num=Vtx_Num();
	const int ele_num=Ele_Num();
		
	Timer timer2;
	timer2.Begin_Loop();
	for(int ele_idx=0; ele_idx <ele_num; ele_idx++){
		MatrixD grad_s; MatrixD hess_s[d][d];
		Stretch_Force(ele_idx, grad_s, hess_s);
		timer2.Record("calculate streching force");
		for(int j=0;j<d;j++){
			Add_Block(b, E()[ele_idx][j], -dt*grad_s.col(j));
			for(int k=0;k<d;k++){
				Add_Block(b,E()[ele_idx][j],-dt*(dt+damping)*hess_s[j][k]*V()[E()[ele_idx][k]]);
				Add_Block_Helper(A, E()[ele_idx][j], E()[ele_idx][k], dt * (dt + damping) * hess_s[j][k]);
			}
		}
		timer2.Record("assemble strething matrix");
	}
	timer2.Output_Profile(std::cout);
	Info("Assemble linear system for stretching", timer.Lap_Time(PhysicalUnits::ms));
}

template<int d> void SoftBodyNonlinearFemThinShell<d>::Update_Implicit_Bending(const real dt) {
	Timer timer;
	const int vtx_num = Vtx_Num();
	const int ele_num = Ele_Num();

	timer.Begin_Loop();
	if constexpr (d == 3) {
		for (int i = 0; i < edges.size(); i++) {
			Eigen::Matrix<real, d, d + 1> grad_b;
			Eigen::Matrix<real, d, d> hess_b[d + 1][d + 1];
			ArrayF<int, d + 1> vtx_idx;
			ArrayF<int, 2> ele_idx;
			if (Junction_Info(i, vtx_idx, ele_idx)) {
				Bend_Force_Approx(i, grad_b, vtx_idx, ele_idx, hess_b);
				timer.Record("calculate bending force");
				for (int j = 0; j < d + 1; j++) {
					Add_Block(b, vtx_idx[j], -dt * grad_b.col(j));

					for (int k = 0; k < d + 1; k++) {
						Add_Block(b, vtx_idx[j], -dt * (dt + damping) * hess_b[j][k] * V()[vtx_idx[k]]);
						Add_Block_Helper(A, vtx_idx[j], vtx_idx[k], dt * (dt + damping) * hess_b[j][k]);
					}
				}
				timer.Record("assemble bending matrix");
			}
		}
	}
	timer.Output_Profile(std::cout);
	Info("Assemble linear system for bending", timer.Lap_Time(PhysicalUnits::ms));
}

template<int d> void SoftBodyNonlinearFemThinShell<d>::Update_Implicit_Boundary_Condition(const real dt){
	Timer timer;
	const int vtx_num=Vtx_Num();
	const int ele_num=Ele_Num();

	for (auto& bc_d : bc.psi_D_values) {
		int node = bc_d.first; VectorD dis = bc_d.second;
		for (int axis = 0; axis < d; axis++) {
			int idx = node * d + axis;
			NonlinearFemFunc<d>::Set_Dirichlet_Boundary_Helper(A, b, idx, dis[axis]);
		}
	}
}

template<int d> void SoftBodyNonlinearFemThinShell<d>::Grad_Stretch(real area, const MatrixD& stress, const MatrixD& x_hat, const MatrixD& Dm_inv, MatrixD& grad) {
	MatrixD P=area*stress;
	MatrixD C_x=C_c()*Dm_inv;
	grad=x_hat*C_x*P*C_x.transpose();
}

template<int d> void SoftBodyNonlinearFemThinShell<d>::Grad_Stretch(const int ele_idx, MatrixD& grad_s) {
	ElasticParam& material = materials[material_id[ele_idx]];
	real avg_h = 0; //average thickness on an element, not considering weighing yet

	ArrayF<VectorD, d> vtx;
	MatrixD x_hat;
	for (int j = 0; j < d; j++) {
		vtx[j] = X()[E()[ele_idx][j]];
		x_hat.col(j) = vtx[j];
		avg_h += hs[E()[ele_idx][j]];
	}

	avg_h *= (real)1 / (real)d;
	real ks = Ks(material.youngs_modulus, avg_h, material.poisson_ratio);

	MatrixD ds;
	NonlinearFemFunc<d>::D(vtx, ds);
	MatrixD deformation = ds * Dm_inv[ele_idx];
	MatrixD strain;
	Deformation_To_Strain(deformation, strain);
	MatrixD stress;
	Strain_To_Stress(ks, material.poisson_ratio, strain, stress);

	Grad_Stretch(areas_hat[ele_idx], stress, x_hat, Dm_inv[ele_idx], grad_s);
}


template<int d> void SoftBodyNonlinearFemThinShell<d>::Stretch_Force(const real& area,const real& ks, const real& poisson_ratio, const MatrixD& stress, const MatrixD& x_hat, const MatrixD& Dm_inv, MatrixD& grad, MatrixD hess[d][d]) {
	MatrixD P = area * stress;
	MatrixD C_x = C_c() * Dm_inv;
	grad = x_hat * C_x * P * C_x.transpose();

	MatrixD R = C_x * P * C_x.transpose();
	MatrixD Q = C_x * C_x.transpose();
	MatrixD y_hat = x_hat * Q;
	MatrixD z_hat = y_hat * x_hat.transpose();

	for (int i = 0; i < d; i++) {
		for (int j = 0; j < d; j++) {
			for (int p = 0; p < d; p++) {
				for (int q = 0; q < d; q++) {
					hess[j][q](i, p) = ((i == p) ? R(q, j) : (real)0) + area * ks * poisson_ratio * y_hat(i, j) * y_hat(p, q) + (real)0.5 * area * ks * (1 - poisson_ratio) * (z_hat(p, i) * Q(q, j) + y_hat(p, j) * y_hat(i, q));
				}
			}
		}
	}
}

template<int d> void SoftBodyNonlinearFemThinShell<d>::Stretch_Force(const int ele_idx, MatrixD& grad_s, MatrixD hess_s[d][d]) {
	ElasticParam& material = materials[material_id[ele_idx]];
	real avg_h = 0; //average thickness on an element, not considering weighing yet

	ArrayF<VectorD, d> vtx;
	MatrixD x_hat;
	for (int j = 0; j < d; j++) {
		vtx[j] = X()[E()[ele_idx][j]];
		x_hat.col(j) = vtx[j];
		avg_h += hs[E()[ele_idx][j]];
	}

	avg_h *= (real)1 / (real)d;
	real ks = Ks(material.youngs_modulus, avg_h, material.poisson_ratio);

	MatrixD ds;
	NonlinearFemFunc<d>::D(vtx, ds);
	MatrixD deformation = ds * Dm_inv[ele_idx];
	MatrixD strain;
	Deformation_To_Strain(deformation, strain);
	MatrixD stress;
	Strain_To_Stress(ks, material.poisson_ratio, strain, stress);

	Stretch_Force(areas_hat[ele_idx], ks, material.poisson_ratio, stress, x_hat, Dm_inv[ele_idx], grad_s, hess_s);
}

template<int d> void SoftBodyNonlinearFemThinShell<d>::Grad_Bend(const Eigen::Matrix<real,d,d+1>& dtheta_dx, real theta, real theta_hat, real lambda, Eigen::Matrix<real,d,d+1>& dE_dx){
	dE_dx=(real)2*lambda*(theta-theta_hat)*dtheta_dx;
}

template<int d> void SoftBodyNonlinearFemThinShell<d>::Grad_Bend(int jct_idx, Eigen::Matrix<real,d,d+1>& grad, const ArrayF<int,d+1>& vtx_idx, const ArrayF<int, 2>& ele_idx) {
	//no implementation
}

template<> void SoftBodyNonlinearFemThinShell<2>::Grad_Bend(int jct_idx, Eigen::Matrix<real,2,3>& grad, const ArrayF<int,3>& vtx_idx, const ArrayF<int, 2>& ele_idx) {
	return; //not implemented yet
}

template<> void SoftBodyNonlinearFemThinShell<3>::Grad_Bend(int jct_idx, Eigen::Matrix<real,3,4>& grad, const ArrayF<int,4>& vtx_idx, const ArrayF<int, 2>& ele_idx) {
	Vector3 n0=Triangle<3>::Normal(X()[E()[ele_idx[0]][0]],X()[E()[ele_idx[0]][1]],X()[E()[ele_idx[0]][2]]);
	Vector3 n1=Triangle<3>::Normal(X()[E()[ele_idx[1]][0]],X()[E()[ele_idx[1]][1]],X()[E()[ele_idx[1]][2]]);

	real theta= Dihedral_Angle(n0,n1,X()[vtx_idx[0]],X()[vtx_idx[1]]);

	Vector2 w0=Barycentric_Weights(X()[vtx_idx[2]],X()[vtx_idx[0]],X()[vtx_idx[1]]);
	Vector2 w1=Barycentric_Weights(X()[vtx_idx[3]],X()[vtx_idx[0]],X()[vtx_idx[1]]);
	real h0=Distance(X()[vtx_idx[2]],X()[vtx_idx[0]],X()[vtx_idx[1]]);
	real h1=Distance(X()[vtx_idx[3]],X()[vtx_idx[0]],X()[vtx_idx[1]]);

	Eigen::Matrix<real,3,4> dtheta;
	dtheta.col(0)=-(w0[0]*n0/h0+w1[0]*n1/h1);
	dtheta.col(1)=-(w0[1]*n0/h0+w1[1]*n1/h1);
	dtheta.col(2)=n0/h0;
	dtheta.col(3)=n1/h1;
	grad=(real)2*lambdas[jct_idx]*(theta-theta_hats[jct_idx])*dtheta;
}

template<> void SoftBodyNonlinearFemThinShell<2>::Grad_Theta(Eigen::Matrix<real, 2, 3>& grad_theta, const ArrayF<int, 3>& vtx_idx, const ArrayF<int, 2>& ele_idx) {
	return; //no implementation
}

template<> void SoftBodyNonlinearFemThinShell<3>::Grad_Theta(Eigen::Matrix<real, 3, 4>& grad_theta, const ArrayF<int, 4>& vtx_idx, const ArrayF<int, 2>& ele_idx) {
	Vector3 n0 = Triangle<3>::Normal(X()[E()[ele_idx[0]][0]], X()[E()[ele_idx[0]][1]], X()[E()[ele_idx[0]][2]]);
	Vector3 n1 = Triangle<3>::Normal(X()[E()[ele_idx[1]][0]], X()[E()[ele_idx[1]][1]], X()[E()[ele_idx[1]][2]]);

	real theta = Dihedral_Angle(n0, n1, X()[vtx_idx[0]], X()[vtx_idx[1]]);

	Vector2 w0 = Barycentric_Weights(X()[vtx_idx[2]], X()[vtx_idx[0]], X()[vtx_idx[1]]);
	Vector2 w1 = Barycentric_Weights(X()[vtx_idx[3]], X()[vtx_idx[0]], X()[vtx_idx[1]]);
	real h0 = Distance(X()[vtx_idx[2]], X()[vtx_idx[0]], X()[vtx_idx[1]]);
	real h1 = Distance(X()[vtx_idx[3]], X()[vtx_idx[0]], X()[vtx_idx[1]]);

	grad_theta.col(0) = -(w0[0] * n0 / h0 + w1[0] * n1 / h1);
	grad_theta.col(1) = -(w0[1] * n0 / h0 + w1[1] * n1 / h1);
	grad_theta.col(2) = n0 / h0;
	grad_theta.col(3) = n1 / h1;
}

template<> void SoftBodyNonlinearFemThinShell<2>::Bend_Force_Approx(int edge_idx, Eigen::Matrix<real,2,3>& grad, const ArrayF<int,3>& vtx_idx, const ArrayF<int, 2>& ele_idx, MatrixD hess[3][3]) {
	return; //not implemented yet
}

template<> void SoftBodyNonlinearFemThinShell<3>::Bend_Force_Approx(int edge_idx, Eigen::Matrix<real,3,4>& grad, const ArrayF<int,4>& vtx_idx, const ArrayF<int, 2>& ele_idx, MatrixD hess[4][4]) {
	Vector3 n0=Triangle<3>::Normal(X()[E()[ele_idx[0]][0]],X()[E()[ele_idx[0]][1]],X()[E()[ele_idx[0]][2]]);
	Vector3 n1=Triangle<3>::Normal(X()[E()[ele_idx[1]][0]],X()[E()[ele_idx[1]][1]],X()[E()[ele_idx[1]][2]]);

	real theta=Dihedral_Angle(n0,n1,X()[vtx_idx[0]],X()[vtx_idx[1]]);

	Vector2 w0=Barycentric_Weights(X()[vtx_idx[2]],X()[vtx_idx[0]],X()[vtx_idx[1]]);
	Vector2 w1=Barycentric_Weights(X()[vtx_idx[3]],X()[vtx_idx[0]],X()[vtx_idx[1]]);
	real h0=Distance(X()[vtx_idx[2]],X()[vtx_idx[0]],X()[vtx_idx[1]]);
	real h1=Distance(X()[vtx_idx[3]],X()[vtx_idx[0]],X()[vtx_idx[1]]);

	Eigen::Matrix<real,3,4> dtheta;
	dtheta.col(0)=-(w0[0]*n0/h0+w1[0]*n1/h1);
	dtheta.col(1)=-(w0[1]*n0/h0+w1[1]*n1/h1);
	dtheta.col(2)=n0/h0;
	dtheta.col(3)=n1/h1;

	grad=(real)2*lambdas[edge_idx]*(theta-theta_hats[edge_idx])*dtheta;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			hess[i][j]=(real)2*lambdas[edge_idx]*(dtheta.col(i))*(dtheta.col(j).transpose());
		}
	}
}

template<> void SoftBodyNonlinearFemThinShell<2>::Bend_Force(int jct_idx, Eigen::Matrix<real, 2, 3>& grad, const ArrayF<int, 3>& vtx_idx, const ArrayF<int, 2>& ele_idx, MatrixD hess[3][3]) {
	return; //not implemented yet
}

template<> void SoftBodyNonlinearFemThinShell<3>::Bend_Force(int jct_idx, Eigen::Matrix<real, 3, 4>& grad, const ArrayF<int, 4>& vtx_idx, const ArrayF<int, 2>& ele_idx, MatrixD hess[4][4]) {
	Vector3 n0 = Triangle<3>::Normal(X()[E()[ele_idx[0]][0]], X()[E()[ele_idx[0]][1]], X()[E()[ele_idx[0]][2]]);
	Vector3 n1 = Triangle<3>::Normal(X()[E()[ele_idx[1]][0]], X()[E()[ele_idx[1]][1]], X()[E()[ele_idx[1]][2]]);

	real theta = Dihedral_Angle(n0, n1, X()[vtx_idx[0]], X()[vtx_idx[1]]);

	ArrayF<int, 4> map0, map1; //vtx_idx to triangle index
	map0.fill((int)-1); map1.fill((int)-1);
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 4; j++) {
			if (E()[ele_idx[0]][i] == vtx_idx[j]) {
				map0[j] = i;
			}

			if (E()[ele_idx[1]][i] == vtx_idx[j]) {
				map1[j] = i;
			}
		}
	}

	real ht0 = hts[ele_idx[0]][map0[2]];
	real ht1 = hts[ele_idx[1]][map1[3]];

	Eigen::Matrix<real, 3, 4> grad_theta;
	const Vector2& r0 = rs[ele_idx[0]][map0[2]];
	const Vector2& r1 = rs[ele_idx[1]][map1[3]];
	grad_theta.col(0) = -(r0[0]*n0  + r1[1]*n1); // reverse order for the second triangle because the shared points should be iterated reversely
	grad_theta.col(1) = -(r0[1]*n0  + r1[0]*n1);
	grad_theta.col(2) = n0 / ht0;
	grad_theta.col(3) = n1 / ht1;

	grad = (real)2 * lambdas[jct_idx] * (theta - theta_hats[jct_idx]) * grad_theta;

	Matrix3 hess_theta[4][4];
	hess_theta[2][3] = hess_theta[3][2] = Matrix3::Zero();

	ArrayF<Vector3, 3>vtx0 {X()[E()[ele_idx[0]][0]], X()[E()[ele_idx[0]][1]], X()[E()[ele_idx[0]][2]] };
	ArrayF<Vector3, 3>vtx1{ X()[E()[ele_idx[1]][0]], X()[E()[ele_idx[1]][1]], X()[E()[ele_idx[1]][2]] };
	real q0 = (real)1 / hts[ele_idx[0]][map0[2]];
	real q1 = (real)1 / hts[ele_idx[1]][map1[3]];

	ArrayF<Matrix3, 3>& grad_n0 = grad_ns[ele_idx[0]];
	ArrayF<Matrix3, 3>& grad_n1 = grad_ns[ele_idx[1]];
	
	int i0_p = map0[2]; //mapped i in triangle index
	for (int j = 0; j <= 2; j++) {
		int j_p = map0[j]; //mapped j
		Vector3 grad_q0;
		Grad_Q(vtx0, i0_p, j_p, ps[ele_idx[0]], q0, ls[ele_idx[0]], areas[ele_idx[0]], grad_q0);
		hess_theta[j][2] = grad_q0 * n0.transpose() + q0 * grad_n0[j_p];
		if (2 != j) {hess_theta[2][j] = hess_theta[j][2].transpose();}
	}

	int i1_p = map1[3];
	for (int j = 0; j <= 3; j++) {
		if (j == 2) { continue; }
		int j_p = map1[j];
		Vector3 grad_q1;
		Grad_Q(vtx1, i1_p, j_p, ps[ele_idx[1]], q1, ls[ele_idx[1]], areas[ele_idx[1]], grad_q1);
		hess_theta[j][3] = grad_q1 * n1.transpose() + q1 * grad_n1[j_p];
		if (3 != j) {hess_theta[3][j] = hess_theta[j][3].transpose();}
	}

	
	ArrayF<ArrayF<Vector3, 2>, 3>& grad_r0 = grad_rs[ele_idx[0]][i0_p];
	ArrayF<ArrayF<Vector3, 2>, 3>& grad_r1 = grad_rs[ele_idx[1]][i1_p];

	hess_theta[0][0] = -grad_r0[map0[0]][0] * n0.transpose() - r0[0] * grad_n0[map0[0]] - grad_r1[map1[0]][1] * n1.transpose() - r1[1] * grad_n1[map1[0]];
	hess_theta[1][1] = -grad_r0[map0[1]][1] * n0.transpose() - r0[1] * grad_n0[map0[1]] - grad_r1[map1[1]][0] * n1.transpose() - r1[0] * grad_n1[map1[1]];
	
	hess_theta[1][0] = -grad_r0[map0[1]][0] * n0.transpose() - r0[0] * grad_n0[map0[1]] - grad_r1[map1[1]][1] * n1.transpose() - r1[1] * grad_n1[map1[1]];
	hess_theta[0][1] = hess_theta[1][0].transpose();

	/*std::cout << "hess_theta[" << jct_idx << "]" << std::endl;
	Array2DF<Matrix3, 4, 4> hess_theta_n = Numerical_Hess_Theta(vtx_idx, ele_idx, grad_theta);
	for (int i = 0; i < 4;i++) {
		for (int j = 0; j < 4;j++) {
			if (!hess_theta[i][j].isApprox(hess_theta_n[i][j], (real)1e-2)) {
				std::cout << "hess_theta[" << i << "][" << j << "] \n" << hess_theta[i][j] << "\n" << "hess_theta_n[" << i << "][" << j << "] \n" << hess_theta_n[i][j] << std::endl;
			}
		}
	}*/

	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			hess[i][j] = (real)2 * lambdas[jct_idx] * ((theta - theta_hats[jct_idx]) * hess_theta[i][j] + (grad_theta.col(i)) * (grad_theta.col(j).transpose()));
		}
	}
}

template<int d> void SoftBodyNonlinearFemThinShell<d>::Update_Bending_Hess_Variables() {
	if constexpr (d == 3) {
#pragma omp parallel for
		for (int ele_idx = 0; ele_idx < Ele_Num(); ele_idx++) {
			ArrayF<VectorD, 3> vtx;
			for (int i = 0; i < 3; i++) { vtx[i] = X()[E()[ele_idx][i]]; }
			for (int i = 0; i < 3; i++) { ls[ele_idx][i] = vtx[(i + 1) % 3] - vtx[(i + 2) % 3]; }
			for (int i = 0; i < 3; i++) { l_norms[ele_idx][i] = ls[ele_idx][i].norm(); }
			for (int i = 0; i < 3; i++) { areas[ele_idx] = Triangle<3>::Area(vtx[0], vtx[1], vtx[2]); }

			for (int i = 0; i < 3; i++) {
				ps[ele_idx][i] = ls[ele_idx][(i + 1) % 3].dot(ls[ele_idx][(i + 1) % 3]) + ls[ele_idx][(i + 2) % 3].dot(ls[ele_idx][(i + 2) % 3]) - ls[ele_idx][i].dot(ls[ele_idx][i]);
			}

			for (int i = 0; i < 3; i++) { hts[ele_idx][i] = (real)2 * areas[ele_idx] / l_norms[ele_idx][i]; }

			for (int i = 0; i < 3; i++) {
				rs[ele_idx][i][0] = ps[ele_idx][(i + 2) % 3] / (real)4 / areas[ele_idx] / l_norms[ele_idx][i];
				rs[ele_idx][i][1] = ps[ele_idx][(i + 1) % 3] / (real)4 / areas[ele_idx] / l_norms[ele_idx][i];
				//method2
				/*Vector2 w = Barycentric_Weights(vtx[i], vtx[(i + 1) % 3], vtx[(i + 2) % 3]);
				rs[ele_idx][i][0] = w[0] / hts[ele_idx][i];
				rs[ele_idx][i][1] = w[1] / hts[ele_idx][i];*/
			}

			Grad_N(vtx, ls[ele_idx], areas[ele_idx], grad_ns[ele_idx]);
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++) {
					Grad_R(vtx, i,j,ps[ele_idx], ls[ele_idx], areas[ele_idx], grad_rs[ele_idx][i][j]);
				}
			}
		}
	}
}


template<int d> inline void SoftBodyNonlinearFemThinShell<d>::Reset_To_Rest_Position() {
#pragma omp parallel for
	for (int i = 0; i < Vtx_Num(); i++) {X()[i] = X0[i];}
}

//Helper functions
template<int d> inline void SoftBodyNonlinearFemThinShell<d>::Strain_To_Stress(real ks, real poisson_ratio, const MatrixD& strain,MatrixD& stress)
{
	stress=ks*((1-poisson_ratio)*strain+poisson_ratio*strain.trace()*MatrixD::Identity());
}

template<int d> inline void SoftBodyNonlinearFemThinShell<d>::Deformation_To_Strain(const MatrixD& deformation, MatrixD& strain) {
	strain=(deformation.transpose()*deformation - MatrixD::Identity()) * (real)0.5;
}

template<int d> const Matrix<real,d> SoftBodyNonlinearFemThinShell<d>::C_c() const{/*not impl*/return MatrixD::Identity();}

template<> const Matrix2 SoftBodyNonlinearFemThinShell<2>::C_c() const{
	Matrix2 C_c;
	C_c<<Vector2((real)-1,(real)1),Vector2::Zero();
	return C_c;
}

template<> const Matrix3 SoftBodyNonlinearFemThinShell<3>::C_c() const{
	Matrix3 C_c;
	C_c<<Vector3((real)-1,(real)1,(real)0),Vector3((real)-1,(real)0,(real)1),Vector3::Zero();
	return C_c;
}

template<int d> inline real SoftBodyNonlinearFemThinShell<d>::Stretching_Energy(real a, real ks, real nu, const Matrix<real,d>& strain) {
	real c= (strain.transpose() * strain).trace();
	real b = pow(strain.trace(), 2);
	return (real)0.5*a*ks*(((real)1-nu)*(strain.transpose()*strain).trace() + nu*pow(strain.trace(),2));
}

template<int d> inline real SoftBodyNonlinearFemThinShell<d>::Stretching_Energy(int ele_idx) {
	ElasticParam& material = materials[material_id[ele_idx]];
	real avg_h = 0; //average thickness on an element, not considering weighing yet

	ArrayF<VectorD, d> vtx;
	MatrixD x_hat;
	for (int j = 0; j < d; j++) {
		vtx[j] = X()[E()[ele_idx][j]];
		x_hat.col(j) = vtx[j];
		avg_h += hs[E()[ele_idx][j]];
	}

	avg_h *= (real)1 / (real)d;
	real ks = Ks(material.youngs_modulus, avg_h, material.poisson_ratio);

	MatrixD ds;
	NonlinearFemFunc<d>::D(vtx, ds);
	MatrixD deformation = ds * Dm_inv[ele_idx];
	MatrixD strain;
	Deformation_To_Strain(deformation, strain);

	return Stretching_Energy(areas[ele_idx],ks, material.poisson_ratio,strain);
}

template<int d> inline real SoftBodyNonlinearFemThinShell<d>::Bending_Energy(real lambda, real theta, real theta_hat) {
	return lambda*(theta-theta_hat)*(theta-theta_hat);
}

template<int d> inline real SoftBodyNonlinearFemThinShell<d>::Bending_Energy(int jct_idx, const ArrayF<int, d + 1>& vtx_idx, const ArrayF<int, 2>& ele_idx) {
	if constexpr (d == 3) {
		Vector3 n0 = Triangle<3>::Normal(X()[E()[ele_idx[0]][0]], X()[E()[ele_idx[0]][1]], X()[E()[ele_idx[0]][2]]);
		Vector3 n1 = Triangle<3>::Normal(X()[E()[ele_idx[1]][0]], X()[E()[ele_idx[1]][1]], X()[E()[ele_idx[1]][2]]);
		real theta = ThinShellAuxFunc::Dihedral_Angle(n0, n1, X()[vtx_idx[0]], X()[vtx_idx[1]]);
		return Bending_Energy(lambdas[jct_idx], theta, theta_hats[jct_idx]);
	}
	else {
		return 0; //to be implemented
	}
}

template<int d> real SoftBodyNonlinearFemThinShell<d>::Total_Stretching_Energy() {
	real stretching_energy = 0;
#pragma omp parallel for reduction(+:stretching_energy)
	for (int i = 0; i < Ele_Num(); i++) {
		ElasticParam& material = materials[material_id[i]];
		real avg_h = 0;

		ArrayF<VectorD, d> vtx;
		MatrixD x_hat;
		for (int j = 0; j < d; j++) {
			vtx[j] = X()[E()[i][j]];
			x_hat.col(j) = vtx[j];
			avg_h += hs[E()[i][j]];
		}
		avg_h *= (real)1 / (real)d;
		real ks = Ks(material.youngs_modulus, avg_h, material.poisson_ratio);

		MatrixD ds;
		NonlinearFemFunc<d>::D(vtx, ds);
		MatrixD deformation = ds * Dm_inv[i];
		MatrixD strain;
		Deformation_To_Strain(deformation, strain);

		stretching_energy += Stretching_Energy(areas_hat[i], ks, material.poisson_ratio, strain);
	}
	return stretching_energy;
}

template<int d> real SoftBodyNonlinearFemThinShell<d>::Total_Bending_Energy() {/*No implementation*/ }

template<> real SoftBodyNonlinearFemThinShell<2>::Total_Bending_Energy() { return 0; /*To be implemented*/ }

template<> real SoftBodyNonlinearFemThinShell<3>::Total_Bending_Energy() {
	real bending_energy = 0;
#pragma omp parallel for reduction(+:bending_energy)
	for (int edge_idx = 0; edge_idx < edges.size(); edge_idx++) {
		ArrayF<int, 4> vtx_idx;
		ArrayF<int, 2> ele_idx;
		if (!Junction_Info(edge_idx, vtx_idx, ele_idx)) {continue;}

		Vector3 n0 = Triangle<3>::Normal(X()[E()[ele_idx[0]][0]], X()[E()[ele_idx[0]][1]], X()[E()[ele_idx[0]][2]]);
		Vector3 n1 = Triangle<3>::Normal(X()[E()[ele_idx[1]][0]], X()[E()[ele_idx[1]][1]], X()[E()[ele_idx[1]][2]]);
		real theta = Dihedral_Angle(n0, n1, X()[vtx_idx[0]], X()[vtx_idx[1]]);

		bending_energy+= Bending_Energy(lambdas[edge_idx], theta, theta_hats[edge_idx]);
	}

	return bending_energy;
}

template<int d> real SoftBodyNonlinearFemThinShell<d>::Lambda(const ArrayF<int, d + 1>& vtx_idx, const ArrayF<int, 2>& ele_idx) {
	real l_hat;
	if constexpr (d == 2) { l_hat = (real)1; }
	else { l_hat = (X0[vtx_idx[0]] - X0[vtx_idx[1]]).norm(); } //May be stored

	real a_hat = areas_hat[ele_idx[0]] + areas_hat[ele_idx[1]];
	ElasticParam& material0 = materials[material_id[ele_idx[0]]];
	ElasticParam& material1 = materials[material_id[ele_idx[1]]];
	real avg_h0, avg_h1, avg_h;
	avg_h0 = avg_h1 = avg_h = (real) 0;
	for (int i = 0; i < d + 1; i++) { avg_h += hs[vtx_idx[i]]; } avg_h /= (real)(d + 1);
	for (int i = 0; i < d; i++) { avg_h0 += hs[E()[ele_idx[0]][i]]; } avg_h0 /= (real)d;
	for (int i = 0; i < d; i++) { avg_h1 += hs[E()[ele_idx[1]][i]]; } avg_h1 /= (real)d;
	real ks0 = Ks(material0.youngs_modulus, avg_h0, material0.poisson_ratio);
	real ks1 = Ks(material1.youngs_modulus, avg_h1, material1.poisson_ratio);
	return Lambda(Kb(avg_h, (real)0.5 * (ks0 + ks1)), a_hat, l_hat);
}

template<> inline real SoftBodyNonlinearFemThinShell<2>::Lambda(real kb, real a_hat, real l_hat) {
	return kb / ((real)4 * a_hat); //the other two is divided in the areas_hat adding together
}

template<> inline real SoftBodyNonlinearFemThinShell<3>::Lambda(real kb, real a_hat, real l_hat) {
	return kb * l_hat * l_hat / ((real)4 * a_hat); //the other two is divided in the areas_hat adding together
}

template<int d> Matrix<real,d> SoftBodyNonlinearFemThinShell<d>::Numerical_Grad_Stretch(int ele_idx, const MatrixD& grad_s) {
	MatrixD grad_s_n;
	const real epsilon=1e-10;
	real stretching_energy= Stretching_Energy(ele_idx);

	ArrayF<VectorD, d> vtx;
	for (int j = 0; j < d; j++) {
		vtx[j] = X()[E()[ele_idx][j]];
	}

	for (int col=0; col<d; col++) {
		for (int row = 0; row < d; row++) {
			real p_tmp=vtx[col][row];
			vtx[col][row]+=epsilon;
			real energy_p=Stretching_Energy(ele_idx);
			grad_s_n(row,col)=(energy_p-stretching_energy)/epsilon;
			vtx[col][row]=p_tmp;
		}
	}

	if (!grad_s_n.isApprox(grad_s, (real)1e-2)) {
		std::cout << "grad_s[" << ele_idx << "]: \n" << grad_s << "\n" << "grad_s_n[" << ele_idx << "]: \n" << grad_s_n << std::endl;
	}
	return grad_s_n;
}

template<int d> Array2DF<Matrix<real, d>, d, d> SoftBodyNonlinearFemThinShell<d>::Numerical_Hess_Stretch(int ele_idx, const MatrixD& grad_s, const MatrixD hess_s[d][d]) {
	Array2DF<MatrixD, d, d> hess_s_n;
	const real epsilon = 1e-5;

	ArrayF<VectorD, d> vtx;
	for (int j = 0; j < d; j++) {
		vtx[j] = X()[E()[ele_idx][j]];
	}

	for (int col = 0; col < d; col++) {
		for (int row = 0; row < d; row++) {
			real p_tmp = vtx[col][row];
			vtx[col][row] += epsilon;
			MatrixD grad_cr;
			Grad_Stretch(ele_idx, grad_cr);
			MatrixD hess_n_cr = (grad_cr - grad_s) / epsilon;
			
			//convert the index from hessian between one element to all other elements to
			//hessian between two vectors
			for (int i = 0; i < d; i++) {
				for (int j = 0; j < d; j++) {
					hess_s_n[j][col](i, row) = hess_n_cr(i,j);
				}
			}

			vtx[col][row] = p_tmp;
		}
	}

	for (int r = 0; r < d; r++) {
		for (int c = 0; c < d; c++) {
			if (!hess_s_n[r][c].isApprox(hess_s[r][c], (real)0.05)) {
				std::cout << "hess_s[" << ele_idx << "]" << std::endl;
				std::cout << "hess_s[" << r << "][" << c << "] \n" << hess_s[r][c] << "\n" << "hess_s_n[" << r << "][" << c << "] \n" << hess_s_n[r][c] << std::endl;

				for (int j = 0; j < d; j++) {
					std::cout << "vtx" << j << std::endl;
					std::cout << X0[E()[ele_idx][j]].transpose() << std::endl;
					std::cout << X()[E()[ele_idx][j]].transpose() << std::endl;
				}
			}
		}
	}

	return hess_s_n;
}

template<int d> Eigen::Matrix<real, d, d + 1> SoftBodyNonlinearFemThinShell<d>::Numerical_Grad_Bend(int jct_idx, ArrayF<int, d + 1>& vtx_idx, ArrayF<int, 2>& ele_idx, const Eigen::Matrix<real, d, d + 1>& grad_b) {}

template<> Eigen::Matrix<real, 2, 3> SoftBodyNonlinearFemThinShell<2>::Numerical_Grad_Bend(int jct_idx, ArrayF<int, 3>& vtx_idx, ArrayF<int, 2>& ele_idx, const Eigen::Matrix<real, 2, 3>& grad_b) { /*Not implemented yet*/ Eigen::Matrix<real, 2, 3> grad_n; return grad_n; }

template<> Eigen::Matrix<real,3,4> SoftBodyNonlinearFemThinShell<3>::Numerical_Grad_Bend(int jct_idx, ArrayF<int, 4>& vtx_idx, ArrayF<int, 2>& ele_idx, const Eigen::Matrix<real, 3, 4>& grad_b)
{
	Eigen::Matrix<real,3,4> grad_b_n;
	real epsilon=1e-10;

	Vector3 n0=Triangle<3>::Normal(X()[E()[ele_idx[0]][0]],X()[E()[ele_idx[0]][1]],X()[E()[ele_idx[0]][2]]);
	Vector3 n1=Triangle<3>::Normal(X()[E()[ele_idx[1]][0]],X()[E()[ele_idx[1]][1]],X()[E()[ele_idx[1]][2]]);
	real theta=Dihedral_Angle(n0,n1,X()[vtx_idx[0]],X()[vtx_idx[1]]);
	real bending_energy=Bending_Energy(lambdas[jct_idx], theta, theta_hats[jct_idx]);

	for (int col=0; col<4; col++) {
		for (int row = 0; row < 3; row++) {
			real p_tmp=X()[vtx_idx[col]][row];
			X()[vtx_idx[col]][row]+=epsilon;
			n0=Triangle<3>::Normal(X()[E()[ele_idx[0]][0]],X()[E()[ele_idx[0]][1]],X()[E()[ele_idx[0]][2]]);
			n1=Triangle<3>::Normal(X()[E()[ele_idx[1]][0]],X()[E()[ele_idx[1]][1]],X()[E()[ele_idx[1]][2]]);
			theta=Dihedral_Angle(n0,n1,X()[vtx_idx[0]],X()[vtx_idx[1]]);
			real energy_p=Bending_Energy(lambdas[jct_idx], theta, theta_hats[jct_idx]);
			grad_b_n(row,col)=(energy_p-bending_energy)/epsilon;
			X()[vtx_idx[col]][row]=p_tmp;
		}
	}
	if (!grad_b_n.isApprox(grad_b, (real)1e-2)) {
		std::cout << "grad_b[" << jct_idx << "]: \n" << grad_b << "\n" << "grad_b_n[" << jct_idx << "]: \n" << grad_b_n << std::endl;
	}
	return grad_b_n;
}


template<int d> Array2DF<Matrix<real,d>, d + 1, d + 1> SoftBodyNonlinearFemThinShell<d>::Numerical_Hess_Bend(int jct_idx, const ArrayF<int, d + 1>& vtx_idx, const ArrayF<int, 2>& ele_idx, const Eigen::Matrix<real, d, d + 1>& grad_b, const MatrixD hess_b[d + 1][d + 1]) {
	Array2DF<MatrixD, d + 1, d + 1> hess_b_n;
	const real epsilon = 1e-10;

	for (int col = 0; col < d+1; col++) {
		for (int row = 0; row < d; row++) {
			real p_tmp = X()[vtx_idx[col]][row];
			X()[vtx_idx[col]][row] += epsilon;

			Eigen::Matrix<real, d, d + 1> grad_cr;
			Grad_Bend(jct_idx, grad_cr, vtx_idx, ele_idx);
			
			Eigen::Matrix<real, d, d + 1> hess_n_cr = (grad_cr - grad_b) / epsilon;

			//convert the index from hessian between one element to all other elements to
			//hessian between two vectors
			for (int i = 0; i < d; i++) {
				for (int j = 0; j < d+1; j++) {
					hess_b_n[j][col](i, row) = hess_n_cr(i, j);
				}
			}

			X()[vtx_idx[col]][row] = p_tmp;
		}
	}

	for (int r = 0; r < d + 1; r++) {
		for (int c = 0; c < d + 1; c++) {
			if (!hess_b_n[r][c].isApprox(hess_b[r][c], (real)1e-2)) {
				std::cout << "jct_idx: " << jct_idx << std::endl;
				std::cout << "hess_b[" << r << "][" << c << "] \n" << hess_b[r][c] << "\n" << "hess_b_n[" << r << "][" << c << "] \n" << hess_b_n[r][c] << std::endl;
			}
		}
	}
	return hess_b_n;
}

template<int d> Array2DF<Matrix<real, d>, d + 1, d + 1>  SoftBodyNonlinearFemThinShell<d>::Numerical_Hess_Theta(const ArrayF<int, d + 1>& vtx_idx,const ArrayF<int, 2>& ele_idx, const Eigen::Matrix<real, d, d + 1>& grad_theta) {
	Array2DF<MatrixD, d + 1, d + 1> hess_n;
	const real epsilon = 1e-10;

	for (int col = 0; col < d + 1; col++) {
		for (int row = 0; row < d; row++) {
			real p_tmp = X()[vtx_idx[col]][row];
			X()[vtx_idx[col]][row] += epsilon;

			Eigen::Matrix<real, d, d + 1> grad_cr;
			Grad_Theta(grad_cr, vtx_idx, ele_idx);

			Eigen::Matrix<real, d, d + 1> hess_n_cr = (grad_cr - grad_theta) / epsilon;

			//convert the index from hessian between one element to all other elements to
			//hessian between two vectors
			for (int i = 0; i < d; i++) {
				for (int j = 0; j < d + 1; j++) {
					hess_n[j][col](i, row) = hess_n_cr(i, j);
				}
			}

			X()[vtx_idx[col]][row] = p_tmp;
		}
	}

	return hess_n;
}

template<class T_ARRAY> int Element_Edges(const Vector2i& v,T_ARRAY& edges)
{edges[0][0]=v[0];edges[1][0]=v[1];return 2;}

template<class T_ARRAY> int Element_Edges(const Vector3i& v,T_ARRAY& edges)
{edges[0]=Vector2i(v[0],v[1]);edges[1]=Vector2i(v[1],v[2]);edges[2]=Vector2i(v[2],v[0]);return 3;}



void Grad_Q(const ArrayF<Vector3, 3>& vtx, const int i, const int j, const Vector3& ps, const real& qs_i, const ArrayF<Vector3, 3>& ls, const real& a, Vector3& grad_q) {
	int jr = (j + 1) % 3;
	int jl = (j + 2) % 3;

	real tmp1 = ps[jr] * qs_i / (real)8 / a / a;
	if (i == jr) { tmp1 -= (real)0.5 / (a * ls[i].norm()); }

	real tmp2 = ps[jl] * qs_i / (real)8 / a / a;
	if (i == jl) { tmp2 -= (real)0.5 / (a * ls[i].norm()); }

	grad_q = tmp1 * ls[jr] - tmp2 * ls[jl];
}

void Grad_R(const ArrayF<Vector3, 3>& vtx, const int i, const int j, const Vector3& ps, const ArrayF<Vector3, 3>& ls, const real& a, ArrayF<Vector3, 2>& grad_r) {
	int jr = (j + 1) % 3;
	int jl = (j + 2) % 3;

	Vector3 grad_a = (real)1 / (real)8 / a * (ps[jl] * ls[jl] - ps[jr] * ls[jr]);
	Vector3 grad_l = ((int)(i == jl)) / ls[jl].norm() * ls[jl] - ((int)(i == jr)) / ls[jr].norm() * ls[jr];

	Vector3 coef2 = (grad_a * ls[i].norm() + grad_l * a) / ((real)4 * a * a * ls[i].dot(ls[i]));
	real coef1 = a * ls[i].norm() * (real)2;

	int i_p = (i + 2) % 3;
	Vector3 tmp1 = ((((int)(i_p == jr)) + ((int)(i_p == j)) - ((int)(i_p == jl))) * ls[jl] + (((int)(i_p == jr)) - ((int)(i_p == j)) - ((int)(i_p == jl))) * ls[jr]);
	grad_r[0] = tmp1 / coef1 - ps[i_p] * coef2;

	i_p = (i + 1) % 3;
	tmp1 = ((((int)(i_p == jr)) + ((int)(i_p == j)) - ((int)(i_p == jl))) * ls[jl] + (((int)(i_p == jr)) - ((int)(i_p == j)) - ((int)(i_p == jl))) * ls[jr]);
	grad_r[1] = tmp1 / coef1 - ps[i_p] * coef2;
}

void Grad_N(const ArrayF<Vector3, 3>& vtx, const ArrayF<Vector3, 3>& ls, const real& a, ArrayF<Matrix3, 3>& grad_n) {
	Vector3 n = Triangle<3>::Normal(vtx[0], vtx[1], vtx[2]);
	for (int i = 0; i < 3; i++) {
		grad_n[i] = -n * ((real)0.5 * ls[i].cross(n) / a).transpose();
	}
}

template<int d> bool SoftBodyNonlinearFemThinShell<d>::Junction_Info(int edge_idx, ArrayF<int, d + 1>& vtx_idx, ArrayF<int, 2>& ele_idx)
{
	if constexpr (d == 3) {
		vtx_idx[0] = edges[edge_idx][0];
		vtx_idx[1] = edges[edge_idx][1];
		Array<int> incident_elements;
		Value_Array(edge_element_hashtable, edges[edge_idx], incident_elements);
		if (incident_elements.size() == 2) { //shared edge
			ele_idx[0] = incident_elements[0], ele_idx[1] = incident_elements[1];
			vtx_idx[2] = Third_Vertex(vtx_idx[0], vtx_idx[1], E()[ele_idx[0]]);
			vtx_idx[3] = Third_Vertex(vtx_idx[0], vtx_idx[1], E()[ele_idx[1]]);

			for (int i = 0; i < 3; i++) {
				if (E()[ele_idx[0]][i] == vtx_idx[2]) { vtx_idx[0] = E()[ele_idx[0]][(i + 1) % 3]; vtx_idx[1] = E()[ele_idx[0]][(i + 2) % 3]; }
			}
			return true;
		}
		else {
			return false;
		}
	}
	else if constexpr (d == 2) {
		//need to be implemented
		//vtx_idx[0] = edges[edge_idx][0];
		//vtx_idx[1] = edges[edge_idx][1];
		//Array<int> incident_elements;
		//Value_Array(edge_element_hashtable, edges[edge_idx], incident_elements);
		//if (incident_elements.size() == 2) { //shared edge
		//	ele_idx[0] = incident_elements[0], ele_idx[1] = incident_elements[1];
		//	vtx_idx[2] = Opposite_Vertex(vtx_idx[0], vtx_idx[1], E()[ele_idx[0]]);
		//	vtx_idx[3] = Opposite_Vertex(vtx_idx[0], vtx_idx[1], E()[ele_idx[1]]);

		//	for (int i = 0; i < 3; i++) {
		//		if (E()[ele_idx[0]][i] == vtx_idx[2]) { vtx_idx[0] = E()[ele_idx[0]][(i + 1) % 3]; vtx_idx[1] = E()[ele_idx[0]][(i + 2) % 3]; }
		//	}
		//	return true;
		//}
		//else {
		//	return false;
		//}
		return false;
	}
}

////////////////////////////////////////////////////////////////////////
//omp accelerated functions

inline void Add_Element_Force_To_Vertices(Array<Vector2>& F,const Vector3i& e,const Matrix2& ff)
{
	Vector2& f0=F[e[0]];
	Vector2& f1=F[e[1]];
	Vector2& f2=F[e[2]];
	Vector2 c0=ff.col(0);
	Vector2 c1=ff.col(1);
	Vector2 c2=c0+c1;
#pragma omp atomic
	f0[0]-=c0[0];
#pragma omp atomic
	f0[1]-=c0[1];
#pragma omp atomic
	f1[0]-=c1[0];
#pragma omp atomic
	f1[1]-=c1[1];
#pragma omp atomic
	f2[0]+=c2[0];
#pragma omp atomic
	f2[1]+=c2[1];
}

inline void Add_Element_Force_To_Vertices(Array<Vector3>& F,const Vector4i& e,const Matrix3& ff)
{ 
	Vector3& f0=F[e[0]];
	Vector3& f1=F[e[1]];
	Vector3& f2=F[e[2]];
	Vector3& f3=F[e[3]];
	Vector3 c0=ff.col(0);
	Vector3 c1=ff.col(1);
	Vector3 c2=ff.col(2);
	Vector3 c3=c0+c1+c2;
#pragma omp atomic
	f0[0]-=c0[0];
#pragma omp atomic
	f0[1]-=c0[1];
#pragma omp atomic
	f0[2]-=c0[2];
#pragma omp atomic
	f1[0]-=c1[0];
#pragma omp atomic
	f1[1]-=c1[1];
#pragma omp atomic
	f1[2]-=c1[2];
#pragma omp atomic
	f2[0]-=c2[0];
#pragma omp atomic
	f2[1]-=c2[1];
#pragma omp atomic
	f2[2]-=c2[2];
#pragma omp atomic
	f3[0]+=c3[0];
#pragma omp atomic
	f3[1]+=c3[1];
#pragma omp atomic
	f3[2]+=c3[2];
}

template<int d> void SoftBodyNonlinearFemThinShell<d>::Add_Block_Helper(SparseMatrix<real>& K, const int i, const int j, const MatrixD& Ks)
{
	SparseFunc::Add_Block<d, MatrixD>(K, i, j, Ks);
}

template<int d> void SoftBodyNonlinearFemThinShell<d>::Set_Block(VectorX& b, const int i, const VectorD& bi)
{
	for (int ii = 0; ii < d; ii++)b[i * d + ii] = bi[ii];
}

template<int d> void SoftBodyNonlinearFemThinShell<d>::Add_Block(VectorX& b, const int i, const VectorD& bi)
{
	for (int ii = 0; ii < d; ii++)b[i * d + ii] += bi[ii];
}

template class SoftBodyNonlinearFemThinShell<2>;
template class SoftBodyNonlinearFemThinShell<3>;
