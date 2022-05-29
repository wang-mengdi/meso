#pragma once
#include "NonlinearFemFunc.h"
#include "SparseMatrixMapping.h"
#include "ConjugateGradient.h"
#include "IOHelper.h"
#include "SparseFunc.h"

//thin shell simulator with extra gradient information for optimization
using namespace Meso;
template<int d> class DiscreteShellTopoOpt:public DiscreteShell<d>
{
	Typedef_VectorD(d); Typedef_MatrixD(d); using Base = DiscreteShell<d>;
public:
	void Output(const bf::path base_path, const int frame) {
		std::string vtu_name = fmt::format("vtu{:04d}.vtu", frame);
		bf::path vtu_path = base_path / bf::path(vtu_name);
		DiscreteShellVTKFunc::Output_VTU_T<d, VectorD>(mesh, hs, vtu_path.string());
	}

	// A = hess
	// b = -grad
	// Adx=b
	void Advance_Quasi_Static(std::function <real(real)> H_To_Rho, int power)
	{
		int iter = 0;
		real err = 1;
		const int vtx_num = Vtx_Num();
		const int ele_num = Ele_Num();
		const real alpha = 0.1;				//use line search for deciding alpha?
		const int max_iter = 1000;
		while (err > 1e-3) {
			if (iter == max_iter) {
				Info("max iteration {} reached!", max_iter);
				break;
			}

			SparseFunc::Set_Value(A, (real)0);
			ArrayFunc::Fill(b,0);

			//add external forces
			for (auto& force : bc.forces) {
				Add_Block(b, force.first, force.second);
			}

			//Stretching
			for (int ele_idx = 0; ele_idx < ele_num; ele_idx++) {
				MatrixD grad_s; MatrixD hess_s[d][d];
				Base::Stretch_Force(ele_idx, grad_s, hess_s);

				//Base::Numerical_Grad_Stretch(ele_idx, grad_s);
				//Base::Numerical_Hess_Stretch(ele_idx, grad_s, hess_s);
				//iterate through verteices in the element
				for (int j = 0; j < d; j++) {
					Add_Block(b, E()[ele_idx][j], pow(H_To_Rho(Base::Ele_H(ele_idx)), power)*-grad_s.col(j));
					for (int k = 0; k < d; k++) {
						SparseFunc::Add_Block<d, MatrixD>(A, E()[ele_idx][j], E()[ele_idx][k], pow(H_To_Rho(Base::Ele_H(ele_idx)), power) * hess_s[j][k]);
					}
				}
			}

			//Bending
			if constexpr (d == 3) {
				Update_Bending_Hess_Variables();
				for (int jct_idx = 0; jct_idx < edges.size(); jct_idx++) {
					Eigen::Matrix<real, d, d + 1> grad_b;
					MatrixD hess_b[d + 1][d + 1];
					ArrayF<int, d + 1> vtx_idx;
					ArrayF<int, 2> ele_idx;
					if (Junction_Info(jct_idx, vtx_idx, ele_idx)) {
						Base::Bend_Force(jct_idx, grad_b, vtx_idx, ele_idx, hess_b);

						//Eigen::Matrix<real, d, d + 1> grad_b_n = Base::Numerical_Grad_Bend(i, theta_hats[i], lambdas[i]);

						for (int j = 0; j < d + 1; j++) {
							Add_Block(b, vtx_idx[j], pow(H_To_Rho(Jct_H(jct_idx)), power)*-grad_b.col(j));
							for (int k = 0; k < d + 1; k++) {
								SparseFunc::Add_Block<d, MatrixD>(A, vtx_idx[j], vtx_idx[k], pow(H_To_Rho(Jct_H(jct_idx)), power) * hess_b[j][k]);
							}
						}
					}
				}
			}

			for (auto& bc_d : bc.psi_D_values) {
				int node = bc_d.first; VectorD dis = bc_d.second;
				for (int axis = 0; axis < d; axis++) {
					int idx = node * d + axis;
					if (iter == 0) { NonlinearFemFunc::Set_Dirichlet_Boundary_Helper(A, b, idx, dis[axis]); }
					else { NonlinearFemFunc::Set_Dirichlet_Boundary_Helper(A, b, idx, (real)0); }
				}
			}

			SparseMatrixMapping<real, DataHolder::DEVICE> meso_mat(A);
			SparseDiagonalPreconditioner<real> meso_sparse_diag_pred(meso_mat);
			ConjugateGradient<real> meso_sparse_cg;
			meso_sparse_cg.Init(&meso_mat, &meso_sparse_diag_pred, false, -1, 1e-5);
			ArrayDv<real> dv_x(dx);
			ArrayDv<real> dv_b(b);

			int iters; real relative_error;
			meso_sparse_cg.Solve(dv_x, dv_b, iters, relative_error);
			//Info("Implicit solve {} iters with relative_error {}", iters, relative_error);
			dx = dv_x;

#pragma omp parallel for
			for (int i = 0; i < vtx_num; i++) {
				for (int j = 0; j < d; j++) {
					X()[i][j] += alpha*dx[i * d + j];
				}
			}

			err = ArrayFunc::Norm<real, DataHolder::HOST>(dx) / dx.size();
			Info("Error is {}", err);
			iter++;
		}
		Info("Quasi_static solve finished with {} Newton iterations: ", iter);
	}

	void Update_Lambdas() {
		int junction_size;
		if constexpr (d == 2) { junction_size = Vtx_Num(); }
		else { junction_size = edges.size(); }

		for (int jct_idx = 0; jct_idx < junction_size; jct_idx++) {
			ArrayF<int, d + 1> vtx_idx;
			ArrayF<int, 2> ele_idx;
			if (Junction_Info(jct_idx, vtx_idx, ele_idx)) { //shared edge
				lambdas[jct_idx] = Base::Lambda(vtx_idx, ele_idx);
			}
		}
	}

	void Allocate_DGrad_DRho(SparseMatrix<real>& dgrad_drho) {
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
					r = ThinShellAuxFunc::Third_Vertex(e[0], e[1], E()[face_idx_0]);
					c = ThinShellAuxFunc::Third_Vertex(e[0], e[1], E()[face_idx_1]);
					Assert(r != -1 && c != -1, "index out of range for finding opposite vertex");
					for (int rr = r * d; rr < (r + 1) * d; rr++)for (int cc = c * d; cc < (c + 1) * d; cc++) {
						triplets.push_back(Triplet<real>(rr, cc, (real)0));
						triplets.push_back(Triplet<real>(cc, rr, (real)0));
					}
				}
			}
		}

		dgrad_drho.setFromTriplets(triplets.begin(), triplets.end());
	}

	//similar calculation as the quasi-static advance step, without Newton's iteration
	void Update_W(Array<real>& b_w, Array<real>& w, SparseMatrix<real>& dgrad_drho, bool update_b_w, std::function<real()> Dh_Drho, std::function <real(real)> H_To_Rho, int power) {
		SparseFunc::Set_Value(A, (real)0);
		SparseFunc::Set_Value(dgrad_drho, (real)0);
		ArrayFunc::Fill(w, 0);
		const int vtx_num = Vtx_Num();
		const int ele_num = Ele_Num();

		//Stretching
		for (int ele_idx = 0; ele_idx < ele_num; ele_idx++) {
			MatrixD grad_s; MatrixD hess_s[d][d]; 
			Stretch_Force(ele_idx, grad_s, hess_s);

			/*MatrixD grad_s_n = Numerical_Grad_Stretch(ele_idx,grad_s);*/
			/*Array2DF<Matrix<real, d>, d, d> hess_s_n = Numerical_Hess_Stretch(ele_idx, grad_s,hess_s);*/

			//iterate through vertices in the element
			for (int j = 0; j < d; j++) {
				int vtx_j_idx = E()[ele_idx][j];
				if (update_b_w) { Add_Block(b_w, vtx_j_idx, pow(H_To_Rho(Base::Ele_H(ele_idx)), power) * grad_s.col(j)); }//Only for E_t
				for (int k = 0; k < d; k++) {
					int vtx_k_idx = E()[ele_idx][k];
					SparseFunc::Add_Block<d, MatrixD>(A, vtx_j_idx, vtx_k_idx, pow(H_To_Rho(Base::Ele_H(ele_idx)), power)*hess_s[j][k]);
					
					if (!bc.Is_Psi_D(vtx_j_idx)) {
						MatrixD ds_drho = MatrixD::Zero();
						VectorD ds_dh = grad_s.col(j) / ((real)d * Base::Ele_H(ele_idx));
						
						/*MatrixD ds_dh_n = Numerical_Ds_Dh(ele_idx, vtx_k_idx, grad_s);
						if (!ds_dh_n.col(j).isApprox(ds_dh, (real)1e-3)) {
							Info("element: {}, vtx: {}", ele_idx, vtx_j_idx);
							Info("ds_h: {}", ds_dh);
							Info("ds_h_n: {}", ds_dh_n.col(j).eval());
						}*/
						
						ds_drho.row(0) = pow(H_To_Rho(Base::Ele_H(ele_idx)), power) * ds_dh * Dh_Drho() + power * pow(H_To_Rho(Base::Ele_H(ele_idx)), power - 1) * grad_s.col(j) / (real)d; //Add zero padding
						SparseFunc::Add_Block<d, MatrixD>(dgrad_drho, vtx_k_idx, vtx_j_idx, ds_drho);
					}
				}
			}
		}

		//Bending
		Update_Bending_Hess_Variables();
		for (int jct_idx = 0; jct_idx < Jct_Num(); jct_idx++) {
			Eigen::Matrix<real, d, d + 1> grad_b;
			MatrixD hess_b[d + 1][d + 1];
			ArrayF<int, d + 1> vtx_idx;
			ArrayF<int, 2> ele_idx;
			if (Junction_Info(jct_idx, vtx_idx, ele_idx)) {
				Bend_Force(jct_idx, grad_b, vtx_idx, ele_idx, hess_b);
				/*Eigen::Matrix<real, d, d + 1> grad_b_n = Numerical_Grad_Bend(jct_idx, vtx_idx, ele_idx,grad_b);*/
				/*Array2DF<Matrix<real, d>, d + 1, d + 1> hess_b_n = Numerical_Hess_Bend(jct_idx, vtx_idx, ele_idx, grad_b, hess_b);*/

				for (int j = 0; j < d + 1; j++) {
					//if (update_b_w) { Add_Block(b_w, vtx_idx[j], pow(H_To_Rho(Jct_H(jct_idx)), power) * grad_b.col(j)); }//Only for E_t
					Vector<real,d+1> dlambda_dh;
					DLambda_Dh(vtx_idx, ele_idx, dlambda_dh);

					for (int k = 0; k < d + 1; k++) {
						SparseFunc::Add_Block<d, MatrixD>(A, vtx_idx[j], vtx_idx[k], pow(H_To_Rho(Jct_H(jct_idx)), power) * hess_b[j][k]);
						if (!bc.Is_Psi_D(vtx_idx[j])) {
							MatrixD db_drho = MatrixD::Zero();
							VectorD db_dh = dlambda_dh[k] * grad_b.col(j) / lambdas[jct_idx];
							
							/*Eigen::Matrix<real, d, d + 1> db_dh_n = Numerical_DGradb_Dh(jct_idx, vtx_idx[k], grad_b, vtx_idx, ele_idx);
							if (!db_dh_n.col(j).isApprox(db_dh, (real)1e-3)) {
								Info("jct: {}, vtx: {}, h_vtx: {}", jct_idx, vtx_idx[j], vtx_idx[k]);
								Info("db_h: {}", db_dh);
								Info("db_h_n: {}", db_dh_n.col(j).eval());
							}*/

							db_drho.row(0) = pow(H_To_Rho(Jct_H(jct_idx)), power) * db_dh * Dh_Drho() + power * pow(H_To_Rho(Jct_H(jct_idx)), power - 1) * grad_b.col(j) / (real)(d - 1); //Add zero padding
							SparseFunc::Add_Block<d, MatrixD>(dgrad_drho, vtx_idx[k], vtx_idx[j], db_drho);
						}
					}
				}
			}
		}

		for (auto& bc_d : bc.psi_D_values) {
			int node = bc_d.first; VectorD dis = bc_d.second;
			for (int axis = 0; axis < d; axis++) {
				int idx = node * d + axis;
				NonlinearFemFunc::Set_Dirichlet_Boundary_Helper(A, b_w, idx, (real)0); //Might need some change here
			}
		}

		//Todo: Boundary condition for dgrad_drho? no need since w is solved to be zero at certain positions

		//w is the same as delta x if the compliance is the same as energy
		/*for (int i = 0; i < particles.Size(); i++) {
			VectorD delta_xi=X()[i] - X0[i];
			for (int dim = 0; dim < d; dim++) {
				w[i * d + dim] = delta_xi[dim];
			}
		}*/

		SparseMatrixMapping<real, DataHolder::DEVICE> meso_mat(A);
		SparseDiagonalPreconditioner<real> meso_sparse_diag_pred(meso_mat);
		ConjugateGradient<real> meso_sparse_cg;
		meso_sparse_cg.Init(&meso_mat, &meso_sparse_diag_pred, false, -1, 1e-5);
		ArrayDv<real> dv_w(w);
		ArrayDv<real> dv_b(b_w);

		int iters; real relative_error;
		meso_sparse_cg.Solve(dv_w, dv_b, iters, relative_error);
		Info("Implicit solve {} iters with relative_error {}", iters, relative_error);
		w = dv_w;
	}

	void Update_DEt_Stretch_DRho(Array<real>& det_stretch_drho, std::function<real()> Dh_Drho, std::function <real(real)> H_To_Rho, int power) {
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
			NonlinearFemFunc::D<d>(vtx, ds);
			MatrixD deformation = ds * Dm_inv[i];
			MatrixD strain;
			Deformation_To_Strain(deformation, strain);

			real e_stretch = Stretching_Energy(areas_hat[i], ks, material.poisson_ratio, strain);
			for (int j = 0; j < d; j++) {
				det_stretch_drho[E()[i][j]] += pow(H_To_Rho(avg_h), power) * e_stretch / (real)d / avg_h * Dh_Drho();
				det_stretch_drho[E()[i][j]] += power * pow(H_To_Rho(avg_h), power - 1) * e_stretch / (real)d;
			}
		}
	}

	void Update_DEt_Bend_DRho(Array<real>& det_bend_drho, std::function<real()> Dh_Drho, std::function <real(real)> H_To_Rho, int power) {
		for (int jct_idx = 0; jct_idx < Jct_Num(); jct_idx++) {
			ArrayF<int, d + 1> vtx_idx;
			ArrayF<int, 2> ele_idx;
			if (!Junction_Info(jct_idx, vtx_idx, ele_idx)) { continue; }
			Vector<real,d+1> dlambda_dh;
			DLambda_Dh(vtx_idx, ele_idx, dlambda_dh);

			//Vector<real, d + 1> dlambda_dh_n = Numerical_DLambda_Dh(jct_idx, vtx_idx, ele_idx,dlambda_dh);

			if constexpr (d == 3) {
				Vector3 n0 = Triangle<3>::Normal(X()[E()[ele_idx[0]][0]], X()[E()[ele_idx[0]][1]], X()[E()[ele_idx[0]][2]]);
				Vector3 n1 = Triangle<3>::Normal(X()[E()[ele_idx[1]][0]], X()[E()[ele_idx[1]][1]], X()[E()[ele_idx[1]][2]]);
				real theta = ThinShellAuxFunc::Dihedral_Angle(n0, n1, X()[vtx_idx[0]], X()[vtx_idx[1]]);
				real theta_d2 = (theta - theta_hats[jct_idx]); theta_d2 *= theta_d2;
				real edge_h = Jct_H(jct_idx);
				
				for (int i = 0; i < d + 1; i++) {
					det_bend_drho[vtx_idx[i]] += pow(H_To_Rho(edge_h), power) * dlambda_dh[i] * theta_d2 * Dh_Drho();
					det_bend_drho[vtx_idx[i]] += power * pow(H_To_Rho(edge_h), power - 1) * Bending_Energy(lambdas[jct_idx], theta, theta_hats[jct_idx]) / (real)(d - 1);
					//Numerical_DEb_Dh(jct_idx, vtx_idx[i], vtx_idx, ele_idx, Bending_Energy(lambdas[jct_idx], theta, theta_hats[jct_idx]), dlambda_dh[i] * theta_d2);
				}
			}
		}
	}

	void DLambda_Dh(const ArrayF<int, d + 1>& vtx_idx, const ArrayF<int, 2>& ele_idx, Vector<real,d+1>& dlambda_dh) {
		real l_hat;
		if constexpr (d == 3) { l_hat = (X0[vtx_idx[0]] - X0[vtx_idx[1]]).norm(); }
		else { l_hat = 1; }
		real a_hat = areas_hat[ele_idx[0]] + areas_hat[ele_idx[1]];
		ElasticParam& mat0 = materials[material_id[ele_idx[0]]];
		ElasticParam& mat1 = materials[material_id[ele_idx[1]]];

		real avg_h0, avg_h1, avg_h;
		avg_h0 = avg_h1 = avg_h = 0;
		for (int i = 0; i < d + 1; i++) { avg_h += hs[vtx_idx[i]]; } avg_h /= (real)(d + 1);
		for (int i = 0; i < d; i++) { avg_h0 += hs[E()[ele_idx[0]][i]]; } avg_h0 /= (real)d;
		for (int i = 0; i < d; i++) { avg_h1 += hs[E()[ele_idx[1]][i]]; } avg_h1 /= (real)d;

		real ks0 = Ks(mat0.youngs_modulus, avg_h0, mat0.poisson_ratio);	//this thickness should use the thickness on face instead?
		real ks1 = Ks(mat1.youngs_modulus, avg_h1, mat1.poisson_ratio);	//although it really simplifies the calculation
		real ks_bar = (real)0.5 * (ks0 + ks1);

		real dks_b_drho0 = (real)1 / (real)6 * mat0.youngs_modulus / ((real)1 - mat0.poisson_ratio * mat0.poisson_ratio) * avg_h * avg_h;
		real dks_b_drho1 = (real)1 / (real)6 * mat1.youngs_modulus / ((real)1 - mat1.poisson_ratio * mat1.poisson_ratio) * avg_h * avg_h;
		real coef1 = l_hat * l_hat / (real)48 / a_hat; // another two is in a_hat adding together two area_hat;
		real coef2 = (real)0.5 * ks_bar * avg_h;

		if constexpr (d == 3) {
			dlambda_dh[0] = coef1 * (dks_b_drho0 + dks_b_drho1 + coef2);
			dlambda_dh[1] = dlambda_dh[0];
			dlambda_dh[2] = coef1 * (dks_b_drho0 + coef2);
			dlambda_dh[3] = coef1 * (dks_b_drho1 + coef2);
		}
	}

	Vector<real,d+1> Numerical_DLambda_Dh(int jct_idx, const ArrayF<int, d + 1>& vtx_idx, const ArrayF<int, 2>& ele_idx, const Vector<real, d + 1>& dlambda_dh) {
		const real dh = 1e-6;
		Vector<real, d + 1> dlambda_dh_n;
		for (int i = 0; i < d + 1; i++) {
			real old_h = hs[vtx_idx[i]];
			real lambda = lambdas[jct_idx];
			hs[vtx_idx[i]] += dh;
			real lambda_new = Lambda(vtx_idx, ele_idx);
			dlambda_dh_n[i] = (lambda_new - lambda) / dh;
			hs[vtx_idx[i]] = old_h;
		}

		if (!dlambda_dh_n.isApprox(dlambda_dh, (real)1e-3)) {
			Info("Junction: {}", jct_idx);
			Info("dlambda_drho: {}", dlambda_dh);
			Info("dlambda_drho_n: {}", dlambda_dh_n);
		}
		return dlambda_dh_n;
	}

	MatrixD Numerical_Ds_Dh(int ele_idx, int v_idx, const MatrixD& grad_s) {
		const real dh = 1e-6;
		MatrixD ds_dh_n;
		real old_h = hs[v_idx];
		hs[v_idx] += dh;
		MatrixD grad_s_new;
		Grad_Stretch(ele_idx, grad_s_new);
		ds_dh_n = (grad_s_new - grad_s) / dh;
		hs[v_idx] = old_h;
		return ds_dh_n;
	}

	Eigen::Matrix<real, d, d + 1> Numerical_DGradb_Dh(int jct_idx, int v_idx, const Eigen::Matrix<real, d, d + 1>& grad_b, const ArrayF<int, 4>& vtx_idx, const ArrayF<int, 2>& ele_idx) {
		const real dh = 1e-6;
		Eigen::Matrix<real, d, d + 1> db_dh_n;
		real old_h = hs[v_idx];
		hs[v_idx] += dh;
		real old_lambda = lambdas[jct_idx];
		lambdas[jct_idx]= Lambda(vtx_idx, ele_idx);
		Eigen::Matrix<real, d, d + 1> grad_b_new;
		Grad_Bend(jct_idx, grad_b_new, vtx_idx, ele_idx);
		db_dh_n = (grad_b_new - grad_b) / dh;
		hs[v_idx] = old_h;
		lambdas[jct_idx] = old_lambda;

		if (!db_dh_n.isApprox(db_dh, (real)1e-3)) {
			Info("Junction: {}", jct_idx);
			Info("db_dh: {}", db_dh);
			Info("db_dh_n: {}", db_dh_n);
		}
		return db_dh_n;
	}

	Eigen::Matrix<real, d, d + 1> Numerical_DGradb_Dlambda(int jct_idx, const Eigen::Matrix<real, d, d + 1>& grad_b, const ArrayF<int, 4>& vtx_idx, const ArrayF<int, 2>& ele_idx) {
		const real dlambda = 1e-6;
		Eigen::Matrix<real, d, d + 1> db_dlambda_n;
		real old_lambda = lambdas[jct_idx];
		lambdas[jct_idx] += dlambda;
		Eigen::Matrix<real, d, d + 1> grad_b_new;
		Grad_Bend(jct_idx, grad_b_new, vtx_idx, ele_idx);
		db_dlambda_n = (grad_b_new - grad_b) / dlambda;
		lambdas[jct_idx] = old_lambda;
		return db_dlambda_n;
	}

	real Numerical_DEb_Dh(int jct_idx,int v_idx, const ArrayF<int, 4>& vtx_idx, const ArrayF<int, 2>& ele_idx, const real Eb, const real dEb_dh) {
		const real dh = 1e-6;
		real old_h = hs[v_idx];
		real old_lambda = lambdas[jct_idx];
		hs[v_idx] += dh;
		lambdas[jct_idx] = Lambda(vtx_idx, ele_idx);
		real Eb_new=Bending_Energy(jct_idx, vtx_idx, ele_idx);
		real dEb_dh_n = (Eb_new - Eb) / dh;
		hs[v_idx] = old_h;
		lambdas[jct_idx] = old_lambda;
		if (abs(dEb_dh_n - dEb_dh) / abs(dEb_dh_n) > 1e-3) {
			Info("jct_idx: {}, v_idx: {}", jct_idx, v_idx);
			Info("Eb_dh: {}, Eb_dh_n: {}", dEb_dh, dEb_dh_n);
		}
		return dEb_dh_n;
	}

	real Total_Bending_Obj(std::function <real(real)> H_To_Rho, int power) {
		real bending_energy = 0;
		if constexpr (d == 3) {
	#pragma omp parallel for reduction(+:bending_energy)
			for (int jct_idx = 0; jct_idx < edges.size(); jct_idx++) {
				ArrayF<int, 4> vtx_idx;
				ArrayF<int, 2> ele_idx;
				if (!Junction_Info(jct_idx, vtx_idx, ele_idx)) { continue; }
				bending_energy += pow(H_To_Rho(Jct_H(jct_idx)), power) * Bending_Energy(jct_idx, vtx_idx, ele_idx);
			}
		}
		return bending_energy;
	}

	real Total_Stretching_Obj(std::function <real(real)> H_To_Rho, int power) {
		real stretching_energy = 0;
	#pragma omp parallel for reduction(+:stretching_energy)
		for (int ele_idx = 0; ele_idx < Ele_Num(); ele_idx++) {
			stretching_energy += pow(H_To_Rho(Base::Ele_H(ele_idx)), power) * Stretching_Energy(ele_idx);
		}
		return stretching_energy;
	}
};