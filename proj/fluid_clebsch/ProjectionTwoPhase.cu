//////////////////////////////////////////////////////////////////////////
// Project a vector field to divergence free
// Copyright (c) (2018-),Bo Zhu, Shiying Xiong
// This file is part of CompleX,whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#include "ProjectionTwoPhase.h"

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
	if (solver_mode == SolverType::AUTO) { Auto_Select_Mode(); }

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
template<int d> void ProjectionTwoPhase<d>::Allocate_System()
{
	if(is_A_initialized&&!update_A)return;
	std::function<bool(const int)> valid_cell=[=](const int idx)->bool{return this->Is_Valid_Cell(mac_grid->grid.Cell_Coord(idx));};
	Build_Grid_Cell_Matrix_Bijective_Mapping(mac_grid->grid,valid_cell,grid_to_matrix,matrix_to_grid);

	////Setup A and div_u
	int n=(int)matrix_to_grid.size();
	A.resize(n,n);p.resize(n);p.fill((real)0);div_u.resize(n);div_u.fill((real)0);
}

template<int d> void ProjectionTwoPhase<d>::Update_A()
{
	if(is_A_initialized&&!update_A) { return; }
	Timer timer; timer.Reset();
	Array<TripletT> elements;
	for (auto r = 0; r < matrix_to_grid.size(); r++) {
		const VectorDi& cell = mac_grid->grid.Cell_Coord(matrix_to_grid[r]);
		////off-diagonal elements
		for (int i = 0; i < Grid<d>::Number_Of_Nb_C(); i++) {
			real term = Off_Diag_Term(cell, i);
			VectorDi nb_cell = Grid<d>::Nb_C(cell, i);
			if (Is_Valid_Cell(nb_cell)) {
				int c = grid_to_matrix(nb_cell);
				elements.push_back(TripletT((int)r, (int)c, term));}}

		////diagonal elements
		real dia_coef = (real)0;
		for (int i = 0; i < Grid<d>::Number_Of_Nb_C(); i++) {
			int axis; VectorDi face; MacGrid<d>::Cell_Incident_Face(cell, i, axis, face);
			dia_coef += Diag_Face_Term(axis, face);}
		elements.push_back(TripletT((int)r, (int)r, dia_coef));}

	if (verbose)timer.Elapse_And_Output_And_Reset("Update A elements");
	A.setFromTriplets(elements.begin(), elements.end());
	A.makeCompressed();
	if (verbose)timer.Elapse_And_Output_And_Reset("Assemble A to sp_mtx");
	is_A_initialized=true;
}

template<int d> void ProjectionTwoPhase<d>::Apply_Jump_Condition_To_b()
{
	if (!use_explicit_surface_tension)return;
	real one_over_dx = (real)1 / mac_grid->grid.dx;
	int b_size = (int)matrix_to_grid.size();
	#pragma omp parallel for
	for (int r = 0; r < b_size; r++) {
		const VectorDi& cell = mac_grid->grid.Cell_Coord(matrix_to_grid[r]);
		const VectorD& pos = mac_grid->grid.Center(cell);
		for (int i = 0; i < Grid<d>::Number_Of_Nb_C(); i++) {
			VectorDi nb_cell = Grid<d>::Nb_C(cell, i);
			if (Is_Valid_Cell(cell) && Is_Valid_Cell(nb_cell)) {
				real phi0 = (*levelset).phi(cell); real phi1 = (*levelset).phi(nb_cell);
				if (LevelSet<d>::Interface(phi0, phi1)) {
					real theta = LevelSet<d>::Theta(phi0, phi1);
					VectorD intf_pos = ((real)1 - theta) * mac_grid->grid.Center(cell) + theta * mac_grid->grid.Center(nb_cell);
					int axis; VectorDi face; MacGrid<d>::Cell_Incident_Face(cell, i, axis, face);
					real rho = (*rho_face)(axis, face);
					real p_sign = phi0 < 0 ? (real)1 : (real)-1;
					div_u[r] += p_sign * Jump_Condition(intf_pos) * one_over_dx / rho;
				}
			}
		}
	}
}

 ////need to specify target_vol (externally) and current_vol 
 template<int d> void ProjectionTwoPhase<d>::Apply_Vol_Control_To_b()
 {
	 if (!use_vol_control)return;
	 if (target_vol == (real)-1) { std::cerr << "[Error] ProjectionTwoPhase: target_vol not set" << std::endl; return; }

	 if (calc_current_vol)current_vol = levelset->Total_Volume();

	 real vol_correction = (target_vol - current_vol) / current_vol;
	 //real cell_vol=pow(mac_grid->grid.dx,d);

	 if (verbose) std::cout << "current volume: " << current_vol << std::endl;

	 if (verbose) std::cout << "vol correction: " << vol_correction << std::endl;

	 int b_size = (int)matrix_to_grid.size();
#pragma omp parallel for
	 for (int r = 0; r < b_size; r++) {
		 const VectorDi& cell = mac_grid->grid.Cell_Coord(matrix_to_grid[r]);
		 if (Is_Fluid_Cell(cell)) {
			 //real vol = levelset->Cell_Fraction(cell) * cell_vol;
			 real cell_div = vol_correction;
			 div_u[r] += vol_control_ks * mac_grid->grid.dx * cell_div;
		 }
	 }
 }

template<int d> void ProjectionTwoPhase<d>::Apply_Implicit_Surface_Tension(const real dt)
{
	if (Is_Interface_Face_Index == nullptr || Dirac == nullptr) { std::cerr << "[Error] ProjectionTwoPhase: Is_Interface_Face_Index or Dirac is null" << std::endl; return; }
	macgrid_to_matrix.Resize(mac_grid->grid.cell_counts, -1);
	matrix_to_macgrid.clear();
	Build_MacGrid_Face_Matrix_Bijective_Mapping(*mac_grid, Is_Interface_Face_Index, macgrid_to_matrix, matrix_to_macgrid);
	int n = matrix_to_macgrid.size();

	SparseMatrixT B;
	VectorX u_new;
	VectorX u_old;

	// setup A, x, and b
	B.resize(n, n);
	u_new.resize(n); u_new.fill((real)0);
	u_old.resize(n); u_old.fill((real)0);
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

#ifdef USE_CUDA
	MultiGridCuda::Preconditioned_Conjugate_Gradient<SparseMatrix<real>, real>(&B, &u_old[0], &u_new[0], 3000, (real)1e-5,/*diagonal precond*/true);	////GPU D-PCG
#else
	KrylovSolver::Params params; KrylovSolver::ICPCG(B, u_new, u_old, params);	////CPU IC-PCG
#endif

#pragma omp parallel for
	for (int r = 0; r < n; r++) {
		const int axis = matrix_to_macgrid[r].first;
		const VectorDi& face = mac_grid->face_grids[axis].Node_Coord(matrix_to_macgrid[r].second);
		(*velocity)(axis, face) = u_new[r];
	}
}

template<int d> void ProjectionTwoPhase<d>::Update_b()
{
	int b_size=(int)matrix_to_grid.size();
	#pragma omp parallel for
	for (auto r = 0; r < b_size; r++) {
		const VectorDi& cell = mac_grid->grid.Cell_Coord(matrix_to_grid[r]);
		real div = (real)0;
		for (int axis = 0; axis < d; axis++) {
			div += (*velocity)(axis, cell + VectorDi::Unit(axis)) - (*velocity)(axis, cell);}
		////Attention: use negative div here. We are solving -lap p=-div u
		div_u[r] = -div;}
}

template<int d> void ProjectionTwoPhase<d>::Correction()
{
	for (int axis = 0; axis < d; axis++) {
		int face_num = mac_grid->face_grids[axis].node_counts.prod();
		#pragma omp parallel for
		for (int i = 0; i < face_num; i++) {
			VectorDi face = mac_grid->face_grids[axis].Node_Coord(i);
			(*velocity)(axis, face) += Velocity_Offset(axis, face);
		}
	}
}

template<int d> void ProjectionTwoPhase<d>::Update_Mat_Id()
{
	if (solver_mode != SolverType::MULTIGRID_AUTO)return;

	multigrid_params.use_irregular_domain=grid_to_matrix.Has(-1);
	if(!multigrid_params.use_irregular_domain)return;

	is_irregular_domain = true;

	mat_id.Resize(mac_grid->grid.cell_counts,0);
	#pragma omp parallel for
	for (int i = 0; i < (int)mat_id.array.size(); i++) {
		if (!Is_Valid_Cell(mac_grid->grid.Cell_Coord(i))) { mat_id.array[i] = -1; }}
}

template<int d> void ProjectionTwoPhase<d>::Build()
{
	Timer timer;timer.Reset();
	Allocate_System();
	if(verbose)timer.Elapse_And_Output_And_Reset("Allocate A");
	Update_A();
	if(verbose)timer.Elapse_And_Output_And_Reset("Assemble A");
	Update_b();
	if(verbose)timer.Elapse_And_Output_And_Reset("Assemble b");
	Apply_Jump_Condition_To_b();
	if(verbose)timer.Elapse_And_Output_And_Reset("Apply jump condition to b");
	Apply_Vol_Control_To_b();
	if(verbose)timer.Elapse_And_Output_And_Reset("Apply volumn control to b");
}

template<int d> void ProjectionTwoPhase<d>::Solve_CPX(void)
{
#ifdef USE_CPX
	if (!cpx_inited) AuxFunc::Crash_With_Info("ProjectionTwoPhase<d>::Solve_CPX not inited");

	int cell_num = mac_grid->grid.cell_counts.prod();
	cell_b.Resize(mac_grid->grid.cell_counts);

	#pragma omp parallel for
	for (int idx = 0; idx < cell_num; idx++) {
		const VectorDi& cell = mac_grid->grid.Cell_Coord(idx);
		int r = grid_to_matrix(cell);
		if (r == -1)cell_b(cell) = 0;
		else cell_b(cell) = div_u[r];}
	std::cout << "after build b" << std::endl;

	cpx_poisson.Update_b(cell_b);
	std::cout << "after update b" << std::endl;

	cpx_poisson.Solve_Fast();
	std::cout << "after solve" << std::endl;

#pragma omp parallel for
	for (auto r = 0; r < matrix_to_grid.size(); r++) {
		const VectorDi& cell = mac_grid->grid.Cell_Coord(matrix_to_grid[r]);
		p[r] = cpx_poisson.x(cell);}
	std::cout << "after x->p" << std::endl;
#else
	AuxFunc::Crash_With_Info("Please compile cpx module with solver mode CPX_GPU");
#endif
}

template<int d> void ProjectionTwoPhase<d>::Solve()
{
	if (solver_mode == SolverType::AUTO) { Auto_Select_Mode(); }
	if (solver_mode == SolverType::KRYLOV_CPU) {
		KrylovSolver::Params params;
		params.verbose = verbose;
		KrylovSolver::ICPCG(A, p, div_u, params);	////CPU IC-PCG
	}
	else if (solver_mode == SolverType::MULTIGRID_AUTO) {
		MultiGrid::Params multigrid_params;
		multigrid_params.use_auto_calculated_levels = true;
		multigrid_params.dof_on_cell = true;
		multigrid_params.block_size = 1;
		multigrid_params.use_color = multigrid_params.use_gpu;
		multigrid_params.use_irregular_domain = true;
		multigrid_params.use_gpu = true;
		Update_Mat_Id();
#ifdef USE_CUDA
		if (multigrid_params.use_gpu) {
			GMGPCG_GPU<d>(A, p, div_u, mac_grid->grid.cell_counts, mat_id, multigrid_params,/*verbose*/false);
		}
		else { GMGPCG_CPU<d>(A, p, div_u, mac_grid->grid.cell_counts, mat_id, multigrid_params); }
#else
		GMGPCG_CPU<d>(A, p, div_u, mac_grid->grid.cell_counts, mat_id, multigrid_params);
#endif
	}
	else if (solver_mode == SolverType::CPX_GPU) { Solve_CPX(); }
}

template<int d> void ProjectionTwoPhase<d>::Project()
{
	Timer timer;					timer.Reset();
	if (use_implicit_surface_tension) Apply_Implicit_Surface_Tension(current_dt);
	Build();						if(verbose)timer.Elapse_And_Output_And_Reset("Build");
	Solve();						if(verbose)timer.Elapse_And_Output_And_Reset("Solve");
	Correction();					if(verbose)timer.Elapse_And_Output_And_Reset("Correction");
}

//////////////////////////////////////////////////////////////////////////
////Check functions
template<int d> void ProjectionTwoPhase<d>::Pressure(Field<real,d>& pressure) const
{
	pressure.Resize(mac_grid->grid.cell_counts,(real)0);
	iterate_cell(iter,mac_grid->grid){const VectorDi& cell=iter.Coord();
		int idx=grid_to_matrix(cell);if(idx==-1)continue;
		pressure(cell)=p[idx];}
}

template<int d> void ProjectionTwoPhase<d>::Pressure_Gradient(FaceField<real,d>& grad_p) const
{
	grad_p.Fill((real)0);
	for (int axis = 0; axis < d; axis++) {
		int face_num = mac_grid->face_grids[axis].node_counts.prod();
		#pragma omp parallel for
		for (int i = 0; i < face_num; i++) {
			VectorDi face = mac_grid->face_grids[axis].Node_Coord(i);
			grad_p(axis, face) = -Velocity_Offset(axis, face);}}
}

template<int d> void ProjectionTwoPhase<d>::Divergence(Field<real,d>& div) const
{
	div.Resize(mac_grid->grid.cell_counts,(real)0);
	int b_size=(int)matrix_to_grid.size();
	#pragma omp parallel for
	for(auto r=0;r<b_size;r++){
		const VectorDi& cell=mac_grid->grid.Cell_Coord(matrix_to_grid[r]);
		real divg=(real)0;
		for(int axis=0;axis<d;axis++){divg+=((*velocity)(axis,cell+VectorDi::Unit(axis))-(*velocity)(axis,cell));}
		div(cell)=divg;}
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
	real cell_p[2]; for (int i = 0; i < 2; i++)cell_p[i] = Is_Valid_Cell(cell[i]) ? p[grid_to_matrix(cell[i])] : (real)0;
	real one_over_dx = (real)1 / mac_grid->grid.dx;
	real p_jump = 0.;
	if (use_explicit_surface_tension){
		if (Is_Valid_Cell(cell[0]) && Is_Valid_Cell(cell[1])) {
			real phi0 = (*levelset).phi(cell[0]); real phi1 = (*levelset).phi(cell[1]);
			if (LevelSet<d>::Interface(phi0, phi1)){
				real theta = LevelSet<d>::Theta(phi0, phi1);
				VectorD intf_pos = ((real)1 - theta) * mac_grid->grid.Center(cell[0]) + theta * mac_grid->grid.Center(cell[1]);
				real rho = (*rho_face)(axis, face);
				real p_sign = phi0 < 0 ? (real)1 : (real)-1;
				p_jump = p_sign * Jump_Condition(intf_pos) * one_over_dx;}}}

	return -((cell_p[1] + p_jump) - cell_p[0]) / (*rho_face)(axis, face);
}

template class ProjectionTwoPhase<2>;
template class ProjectionTwoPhase<3>;