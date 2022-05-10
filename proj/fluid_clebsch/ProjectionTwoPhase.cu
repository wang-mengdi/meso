//////////////////////////////////////////////////////////////////////////
// Project a vector field to divergence free
// Copyright (c) (2018-),Bo Zhu, Shiying Xiong
// This file is part of CompleX,whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#include "ProjectionTwoPhase.h"
#include "SparseFunc.h"
#include "KrylovSolver.h"
#include "SPX_Timer.h"
#include "SparseMatrixMapping.h"

//////////////////////////////////////////////////////////////////////////
////Constructor
template<int d> ProjectionTwoPhase<d>::ProjectionTwoPhase(MacGrid<d>* _mac_grid, FaceField<real, d>* _velocity, FaceField<real, d>* _rho_face, LevelSet<d>* _levelset, Field<ushort, d>* _type, BoundaryConditionMacGrid<d>* _bc)
{
	Initialize(_mac_grid, _velocity, _rho_face, _levelset, _type, _bc);
}

template<int d> ProjectionTwoPhase<d>::~ProjectionTwoPhase()
{
	if (mac_grid != nullptr && own_grid) delete mac_grid;
	if (velocity != nullptr && own_velocity) delete velocity;
	if (rho_face != nullptr && own_rho_face) delete rho_face;
	if (levelset != nullptr && own_levelset) delete levelset;
	if (type != nullptr && own_type) delete type;
	if (bc != nullptr && own_bc) delete bc;
}

template<int d> void ProjectionTwoPhase<d>::Initialize(MacGrid<d>* _mac_grid, FaceField<real, d>* _velocity, FaceField<real, d>* _rho_face, LevelSet<d>* _levelset, Field<ushort, d>* _type, BoundaryConditionMacGrid<d>* _bc)
{
	if (mac_grid != nullptr && own_grid) delete mac_grid;
	if (_mac_grid == nullptr) { mac_grid = new MacGrid<d>(); own_grid = true; }
	else { mac_grid = _mac_grid; own_grid = false; }

	const VectorDi& counts = mac_grid->grid.cell_counts;

	if (velocity != nullptr && own_velocity)delete velocity;
	if (_velocity == nullptr) { velocity = new FaceField<real, d>(counts, (real)0); own_velocity = true; }
	else { velocity = _velocity; own_velocity = false; }

	if (rho_face != nullptr && own_rho_face)delete rho_face;
	if (_rho_face == nullptr) { rho_face = new FaceField<real, d>(counts, (real)0); own_rho_face = true; }
	else { rho_face = _rho_face; own_rho_face = false; }

	if(_levelset!=nullptr&&own_levelset)delete levelset;
	if(_levelset==nullptr){levelset=new LevelSet<d>(mac_grid->grid);own_levelset=true;}
	else{levelset=_levelset;own_levelset=false;}

	if (type != nullptr && own_type)delete type;
	if (_type == nullptr) { type = new Field<ushort, d>(counts, (ushort)CellType::Fluid); own_type = true; }
	else { type = _type; own_type = false; }

	if (bc != nullptr && own_bc)delete bc;
	if (_bc == nullptr) { bc = new BoundaryConditionMacGrid<d>(*_mac_grid); own_bc = true; }
	else { bc = _bc; own_bc = false; }
}

//////////////////////////////////////////////////////////////////////////
////Build linear system

template<int d> void ProjectionTwoPhase<d>::Update_A()
{
	//if (is_A_initialized && !update_A)return;
	//std::function<bool(const int)> valid_cell = [=](const int idx)->bool {return this->Is_Valid_Cell(mac_grid->grid.Cell_Coord(idx)); };
	//Build_Grid_Cell_Matrix_Bijective_Mapping(mac_grid->grid, valid_cell, grid_to_matrix, matrix_to_grid);

	//////Setup A and div_u
	//int n = (int)matrix_to_grid.size();
	//A.resize(n, n); p.resize(n); p.fill((real)0); div_u.resize(n); div_u.fill((real)0);

	Meso::Grid<d> meso_grid(mac_grid->grid.cell_counts);
	meso_fixed_host.Init(meso_grid);
	meso_fixed_host.Calc_Cells(
		[&](const VectorDi cell) {
			return !this->Is_Valid_Cell(cell);
		}
	);
	meso_rho_host.Init(meso_grid);
	meso_rho_host.Calc_Faces(
		[&](const int axis, const VectorDi face)->real {
			VectorDi cell[2]; for (int i = 0; i < 2; i++)cell[i] = MacGrid<d>::Face_Incident_Cell(axis, face, i);
			if (!Is_Valid_Cell(cell[0])) std::swap(cell[0], cell[1]);
			if (!Is_Valid_Cell(cell[0])) return (real)0;
			if (bc->Is_Psi_N(axis, face))return 0;
			return 1.0 / (*rho_face)(axis, face);
		}
	);
	meso_poisson.Init(meso_grid, meso_rho_host, meso_fixed_host);
	meso_mg.Init_Poisson(meso_poisson, 2, 2);
	meso_cg.Init(&meso_poisson, &meso_mg, true, -1, 1e-9);
	//meso_cg.Init(&meso_poisson, nullptr, true, -1, 1e-5);
}

template<int d> void ProjectionTwoPhase<d>::Apply_Jump_Condition_To_b()
{
	//if (!use_explicit_surface_tension)return;
	real one_over_dx = (real)1 / mac_grid->grid.dx;
	meso_div_host.Exec_Nodes(
		[&](const VectorDi cell) {
			if (!Is_Valid_Cell(cell)) return;
			VectorD pos = mac_grid->grid.Center(cell);
			for (int i = 0; i < Grid<d>::Number_Of_Nb_C(); i++) {
				VectorDi nb_cell = Grid<d>::Nb_C(cell, i);
				if (Is_Valid_Cell(nb_cell)) {
					real phi0 = (*levelset).phi(cell); 
					real phi1 = (*levelset).phi(nb_cell);
					if (LevelSet<d>::Interface(phi0, phi1)) {
						real theta = LevelSet<d>::Theta(phi0, phi1);
						VectorD intf_pos = ((real)1 - theta) * mac_grid->grid.Center(cell) + theta * mac_grid->grid.Center(nb_cell);
						int axis; VectorDi face; MacGrid<d>::Cell_Incident_Face(cell, i, axis, face);
						real rho = (*rho_face)(axis, face);
						real p_sign = phi0 < 0 ? (real)1 : (real)-1;
						meso_div_host(cell) += p_sign * Jump_Condition(intf_pos) * one_over_dx / rho;
					}
				}
			}
		}
	);

	//if (!use_explicit_surface_tension)return;
	//real one_over_dx = (real)1 / mac_grid->grid.dx;
	//int b_size = (int)matrix_to_grid.size();
	//#pragma omp parallel for
	//for (int r = 0; r < b_size; r++) {
	//	const VectorDi& cell = mac_grid->grid.Cell_Coord(matrix_to_grid[r]);
	//	const VectorD& pos = mac_grid->grid.Center(cell);
	//	for (int i = 0; i < Grid<d>::Number_Of_Nb_C(); i++) {
	//		VectorDi nb_cell = Grid<d>::Nb_C(cell, i);
	//		if (Is_Valid_Cell(cell) && Is_Valid_Cell(nb_cell)) {
	//			real phi0 = (*levelset).phi(cell); real phi1 = (*levelset).phi(nb_cell);
	//			if (LevelSet<d>::Interface(phi0, phi1)) {
	//				real theta = LevelSet<d>::Theta(phi0, phi1);
	//				VectorD intf_pos = ((real)1 - theta) * mac_grid->grid.Center(cell) + theta * mac_grid->grid.Center(nb_cell);
	//				int axis; VectorDi face; MacGrid<d>::Cell_Incident_Face(cell, i, axis, face);
	//				real rho = (*rho_face)(axis, face);
	//				real p_sign = phi0 < 0 ? (real)1 : (real)-1;
	//				div_u[r] += p_sign * Jump_Condition(intf_pos) * one_over_dx / rho;
	//			}
	//		}
	//	}
	//}
}

 ////need to specify target_vol (externally) and current_vol 
 template<int d> void ProjectionTwoPhase<d>::Apply_Vol_Control_To_b()
 {
	 //if (!use_vol_control)return;
	 if (target_vol == (real)-1) { std::cerr << "[Error] ProjectionTwoPhase: target_vol not set" << std::endl; return; }

	 if (calc_current_vol)current_vol = levelset->Total_Volume();

	 real vol_correction = (target_vol - current_vol) / current_vol;
	 //real cell_vol=pow(mac_grid->grid.dx,d);

	 if (verbose) std::cout << "current volume: " << current_vol << std::endl;

	 if (verbose) std::cout << "vol correction: " << vol_correction << std::endl;

	 meso_div_host.Exec_Nodes(
		 [&](const VectorDi cell) {
			 if (Is_Fluid_Cell(cell)) {
				 real cell_div = vol_correction;
				 meso_div_host(cell) += vol_control_ks * mac_grid->grid.dx * cell_div;
			 }
		 }
	 );

//	 int b_size = (int)matrix_to_grid.size();
//#pragma omp parallel for
//	 for (int r = 0; r < b_size; r++) {
//		 const VectorDi& cell = mac_grid->grid.Cell_Coord(matrix_to_grid[r]);
//		 if (Is_Fluid_Cell(cell)) {
//			 //real vol = levelset->Cell_Fraction(cell) * cell_vol;
//			 real cell_div = vol_correction;
//			 div_u[r] += vol_control_ks * mac_grid->grid.dx * cell_div;
//		 }
//	 }
 }

template<int d> void ProjectionTwoPhase<d>::Apply_Implicit_Surface_Tension(const real dt)
{
	if (Is_Interface_Face_Index == nullptr || Dirac == nullptr) { std::cerr << "[Error] ProjectionTwoPhase: Is_Interface_Face_Index or Dirac is null" << std::endl; return; }
	FaceField<int, d> macgrid_to_matrix;
	Array<std::pair<int, int> > matrix_to_macgrid;
	macgrid_to_matrix.Resize(mac_grid->grid.cell_counts, -1);
	matrix_to_macgrid.clear();
	Build_MacGrid_Face_Matrix_Bijective_Mapping(*mac_grid, Is_Interface_Face_Index, macgrid_to_matrix, matrix_to_macgrid);
	int n = matrix_to_macgrid.size();

	SparseMatrixT B;
	//VectorX u_new;
	//VectorX u_old;

	// setup A, x, and b
	B.resize(n, n);
	//u_new.resize(n); u_new.fill((real)0);
	Meso::Array<real> u_old(n); Meso::ArrayFunc::Fill(u_old, 0);
	//u_old.resize(n); u_old.fill((real)0);
	Array<TripletT> elements;

	for (int r = 0; r < n; r++) {
		const int axis = matrix_to_macgrid[r].first;
		const VectorDi& face = mac_grid->face_grids[axis].Node_Coord(matrix_to_macgrid[r].second);
		VectorD pos = mac_grid->Face_Center(axis, face);
		real kappa = (*levelset).Curvature(pos);
		VectorD normal = (*levelset).Normal(pos);
		real dia_coef = (real)1;
		u_old[r] = (*velocity)(axis, face) - sigma * Dirac((*levelset).Phi(pos)) * dt * kappa * normal[axis];	//neglect Jacobi and Hessian
		for (int i = 0; i < Grid<d>::Number_Of_Nb_C(); i++) {
			VectorDi nb_face = Grid<d>::Nb_C(face, i);
			VectorD nb_pos = mac_grid->Face_Center(axis, nb_face);
			if (!mac_grid->Valid_Face(axis, nb_face) || (*bc).Is_Psi_N(axis, nb_face)) continue;
			real a = sigma * Dirac((*levelset).Phi((pos + nb_pos) * (real).5)) * dt * dt / (mac_grid->grid.dx * mac_grid->grid.dx);
			dia_coef += a;
			int c = macgrid_to_matrix(axis, nb_face);
			if (Is_Interface_Face_Index(std::make_pair(axis, mac_grid->face_grids[axis].Node_Index(nb_face)))) {
				elements.push_back(TripletT(r, c, -a));
			}
			else u_old[r] += a * (*velocity)(axis, nb_face);
		}
		elements.push_back(TripletT(r, r, dia_coef));
	}
	B.setFromTriplets(elements.begin(), elements.end()); B.makeCompressed();

	Meso::SparseMatrixMapping<real, Meso::DataHolder::DEVICE> meso_mat(B);
	Meso::SparseDiagonalPreconditioner<real> meso_sparse_diag_pred(meso_mat);
	Meso::ConjugateGradient<real> meso_sparse_cg;
	meso_sparse_cg.Init(&meso_mat, &meso_sparse_diag_pred, false, -1, 1e-5);
	Meso::ArrayDv<real> u_new_dev(meso_mat.XDoF());
	Meso::ArrayDv<real> u_old_dev = u_old;
	int iters; real relative_error;
	meso_sparse_cg.Solve(u_new_dev, u_old_dev, iters, relative_error);
	Meso::Info("implicit surface tension solve {} iters with relative_error {}", iters, relative_error);
	Meso::Array<real> u_new_host = u_new_dev;

	//KrylovSolver::Params params; KrylovSolver::ICPCG(B, u_new, u_old, params);	////CPU IC-PCG

#pragma omp parallel for
	for (int r = 0; r < n; r++) {
		const int axis = matrix_to_macgrid[r].first;
		const VectorDi& face = mac_grid->face_grids[axis].Node_Coord(matrix_to_macgrid[r].second);
		(*velocity)(axis, face) = u_new_host[r];
	}
}

template<int d> void ProjectionTwoPhase<d>::Update_b()
{
	Meso::Grid<d> meso_grid(mac_grid->grid.cell_counts);
	meso_div_host.Init(meso_grid);
	meso_div_host.Calc_Cells(
		[&](const VectorDi cell)->real {
			if (!mac_grid->grid.Valid_Cell(cell)) return 0;
			real div = (real)0;
			for (int axis = 0; axis < d; axis++) {
				div += (*velocity)(axis, cell + VectorDi::Unit(axis)) - (*velocity)(axis, cell);
			}
			//solve -lap p=-div u
			return -div;
		}
	);

	//meso_velocity_host.Init(meso_grid);
	//meso_velocity_host.Calc_Faces(
	//	[&](const int axis, const VectorDi face) {
	//		if (mac_grid->Valid_Face(axis, face)) return (*velocity)(axis, face);
	//		else return 0;
	//	}
	//);
	//meso_velocity_dev = meso_velocity_dev;
	//Meso::Exterior_Derivative(meso_divergence_dev, meso_velocity_dev);

	//int b_size=(int)matrix_to_grid.size();
	//#pragma omp parallel for
	//for (auto r = 0; r < b_size; r++) {
	//	const VectorDi& cell = mac_grid->grid.Cell_Coord(matrix_to_grid[r]);
	//	real div = (real)0;
	//	for (int axis = 0; axis < d; axis++) {
	//		div += (*velocity)(axis, cell + VectorDi::Unit(axis)) - (*velocity)(axis, cell);}
	//	////Attention: use negative div here. We are solving -lap p=-div u
	//	div_u[r] = -div;}
}

template<int d> void ProjectionTwoPhase<d>::Correction()
{
	for (int axis = 0; axis < d; axis++) {
		int face_num = mac_grid->face_grids[axis].node_counts.prod();
		#pragma omp parallel for
		for (int i = 0; i < face_num; i++) {
			VectorDi face = mac_grid->face_grids[axis].Node_Coord(i);
			VectorDi middle = mac_grid->grid.cell_counts / 2;
			real offset = Velocity_Offset(axis, face);
			for (int i = 0; i < d; i++) { if (i != 1) middle[i] = 0; }
			if (axis == 0 && face == middle) { std::cout << "face vel offset: " << offset << std::endl; }
			(*velocity)(axis, face) += offset;
		}
	}
}

template<int d> void ProjectionTwoPhase<d>::Build()
{
	Meso::Vector<int,d> middle = mac_grid->grid.cell_counts / 2;
	for (int i = 0; i < d; i++) { if (i != 1) middle[i] = 0; }
	Timer timer;timer.Reset();
	Update_A();
	if(verbose)timer.Elapse_And_Output_And_Reset("Assemble A");
	Update_b();
	std::cout << "div after Update b: " << meso_div_host(middle) << std::endl;
	if(verbose)timer.Elapse_And_Output_And_Reset("Assemble b");
	if (use_explicit_surface_tension) {
		Apply_Jump_Condition_To_b();
		std::cout << "div after Jump Condition b: " << meso_div_host(middle) << std::endl;
		if (verbose)timer.Elapse_And_Output_And_Reset("Apply jump condition to b");
	}
	if (use_vol_control) {
		Apply_Vol_Control_To_b();
		if (verbose)timer.Elapse_And_Output_And_Reset("Apply volumn control to b");
	}
}

template<int d> void ProjectionTwoPhase<d>::Solve()
{
	//KrylovSolver::Params params;
	//params.verbose = verbose;
	//KrylovSolver::ICPCG(A, p, div_u, params);	////CPU IC-PCG
	//Solve_CPX();
	meso_div_dev = meso_div_host;
	int iter; real relative_error;
	meso_pressure_dev.Init(meso_div_dev.grid);
	meso_cg.Solve(meso_pressure_dev.Data(), meso_div_dev.Data(), iter, relative_error);
	Meso::Info("MESO solved {} iters with relative error {}", iter, relative_error);
	meso_pressure_host = meso_pressure_dev;
}

template<int d> void ProjectionTwoPhase<d>::Project()
{
	Meso::Vector<int, d> middle = mac_grid->grid.cell_counts / 2;
	for (int i = 0; i < d; i++) { if (i != 1) middle[i] = 0; }
	Timer timer;					timer.Reset();
	if (use_implicit_surface_tension) Apply_Implicit_Surface_Tension(current_dt);
	Build();						if(verbose)timer.Elapse_And_Output_And_Reset("Build");

	//Meso::Info("div_u: \n{}", meso_div_host);

	Solve();						if(verbose)timer.Elapse_And_Output_And_Reset("Solve");
	
	//Meso::Info("solved pressure: \n{}", meso_pressure_dev);
	std::cout << "pressure after solve: " << meso_pressure_host(middle) << std::endl;

	//Meso::FieldDv<float, d> poisson_result;
	//poisson_result.Init(meso_poisson.fixed.grid);
	//meso_poisson.Apply(poisson_result.Data(), meso_pressure_dev.Data());
	//poisson_result -= meso_div_dev;
	//Info("poisson residual max {}", poisson_result.Max_Abs());

	Correction();					if(verbose)timer.Elapse_And_Output_And_Reset("Correction");
	std::cout << "vel after correct: " << (*velocity)(0, middle) << std::endl;
	
	//Meso::FaceField<float, d> meso_velocity_host;
	//meso_velocity_host.Init(meso_pressure_dev.grid);
	//meso_velocity_host.Calc_Faces(
	//	[&](const int axis, const VectorDi face)->float {
	//		if (mac_grid->Valid_Face(axis, face)) return (*velocity)(axis, face);
	//		else return 0;
	//	}
	//);
	//Meso::FaceFieldDv<float, d> meso_velocity_dev;
	//meso_velocity_dev = meso_velocity_dev;
	//Meso::Exterior_Derivative(meso_div_dev, meso_velocity_dev);
	//Meso::Info("after projection max abs: {}", meso_div_dev.Max_Abs());
}

//////////////////////////////////////////////////////////////////////////
////Physical interface functions that defines the problem
//see: https://wmdcstdio.com/2021/07/11/projection-matrix-terms/
template<int d> real ProjectionTwoPhase<d>::Off_Diag_Term(const VectorDi& fluid_cell, const int& nbidx) const
{
	VectorDi nb_cell = Grid<d>::Nb_C(fluid_cell, nbidx);
	int axis; VectorDi face; MacGrid<d>::Cell_Incident_Face(fluid_cell, nbidx, axis, face);
	if (bc->Is_Psi_N(axis, face))return 0;
	if (Is_Valid_Cell(nb_cell)){return -1./(*rho_face)(axis, face);}
	return 0;
}

template<int d> real ProjectionTwoPhase<d>::Diag_Face_Term(const int& axis, const VectorDi& face) const
{
	VectorDi cell[2]; for (int i = 0; i < 2; i++)cell[i] = MacGrid<d>::Face_Incident_Cell(axis, face, i);
	if (!Is_Valid_Cell(cell[0])) std::swap(cell[0], cell[1]);
	if (!Is_Valid_Cell(cell[0])) return (real)0;
	if (bc->Is_Psi_N(axis, face)) return (real)0;
	return 1./(*rho_face)(axis, face);
}

template<int d> real ProjectionTwoPhase<d>::Velocity_Offset(const int& axis, const VectorDi& face) const
{
	if (bc->Is_Psi_N(axis, face)) return 0;
	VectorDi cell[2]; for (int i = 0; i < 2; i++)cell[i] = MacGrid<d>::Face_Incident_Cell(axis, face, i);
	//real cell_p[2]; for (int i = 0; i < 2; i++)cell_p[i] = Is_Valid_Cell(cell[i]) ? p[grid_to_matrix(cell[i])] : (real)0;
	real cell_p[2]; for (int i = 0; i < 2; i++)cell_p[i] = Is_Valid_Cell(cell[i]) ? meso_pressure_host(cell[i]) : (real)0;
	real one_over_dx = (real)1 / mac_grid->grid.dx;
	real p_jump = 0.;
	VectorDi middle = mac_grid->grid.cell_counts / 2;
	for (int i = 0; i < d; i++) { if (i != 1) middle[i] = 0; }
	if (use_explicit_surface_tension) {
		if (Is_Valid_Cell(cell[0]) && Is_Valid_Cell(cell[1])) {
			real phi0 = (*levelset).phi(cell[0]); real phi1 = (*levelset).phi(cell[1]);
			if (LevelSet<d>::Interface(phi0, phi1)) {
				real theta = LevelSet<d>::Theta(phi0, phi1);
				VectorD intf_pos = ((real)1 - theta) * mac_grid->grid.Center(cell[0]) + theta * mac_grid->grid.Center(cell[1]);
				real rho = (*rho_face)(axis, face);
				real p_sign = phi0 < 0 ? (real)1 : (real)-1;
				p_jump = p_sign * Jump_Condition(intf_pos) * one_over_dx;
			}
		}
	}
	if (axis == 0 && face == middle) {
		std::cout << "cell p[1]: " << cell_p[0] << std::endl;
		std::cout << "p_jump: " << p_jump << std::endl;
		std::cout << "cell p[1]: " << cell_p[1] << std::endl;
		std::cout << "rho: " << (*rho_face)(axis, face) << std::endl;
	}
	return -((cell_p[1] + p_jump) - cell_p[0]) / (*rho_face)(axis, face);
}

template class ProjectionTwoPhase<2>;
template class ProjectionTwoPhase<3>;