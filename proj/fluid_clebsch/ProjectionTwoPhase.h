//////////////////////////////////////////////////////////////////////////
// Project a vector field to divergence free on a MAC grid
// Copyright (c) (2018-), Bo Zhu, Shiying Xiong
// This file is part of CompleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////
//// This projection solver currently only takes zero Neumann bc. The velocities on the Neumann boundary need to be set to the correct values before the projection.
//// The calculated pressure value is the real pressure value divided by delta_x.
//////////////////////////////////////////////////////////////////////////

#ifndef __ProjectionTwoPhase_h__
#define __ProjectionTwoPhase_h__
#include <functional>
#include "MacGrid.h"
#include "SPX_FaceField.h"
#include "SPX_Field.h"
#include "BoundaryCondition.h"
#include "GmgPcgSolverCPU.h"
#include "TypeFunc.h"
#include "LevelSet.h"

#ifdef USE_CUDA
#include "GmgPcgSolverGPU.h"
#endif

template<int d> class ProjectionTwoPhase
{Typedef_VectorDii(d);
public:
	////data
	MacGrid<d>* mac_grid;
	FaceField<real,d>* velocity=nullptr;
	FaceField<real,d>* rho_face=nullptr;
    LevelSet<d>* levelset=nullptr;
	Field<ushort,d>* type=nullptr;
	BoundaryConditionMacGrid<d>* bc=nullptr;
	bool own_grid=false;
	bool own_velocity=false;
	bool own_rho_face=false;
	bool own_levelset=false;
	bool own_type=false;
	bool own_bc=false;
    ////flags
    bool use_implicit_surface_tension = false;
    bool use_explicit_surface_tension = true;

	////linear solver
	FaceField<int, d> macgrid_to_matrix;
	Array<std::pair<int, int> > matrix_to_macgrid;

	real sigma = (real)1e-2;					////surface tension coefficient for the default pressure jump
	real current_dt = (real)1;					////need to set dt when using the jump condition because dt is absorbed in p when calculating the projection

    ////callback functions
	//std::function<bool(const int)> Is_Fluid_Interior_Cell_Index = nullptr;
	std::function<real(const VectorD&)> Jump_Condition = std::bind(&ProjectionTwoPhase<d>::Pressure_Jump, this, std::placeholders::_1);
	std::function<bool(const std::pair<int, int>)> Is_Interface_Face_Index = std::bind(&ProjectionTwoPhase<d>::Is_Levelset_Interface_Face_Index, this, std::placeholders::_1);
	std::function<real(real)> Dirac = std::bind(&ProjectionTwoPhase<d>::Levelset_Dirac, this, std::placeholders::_1);;

	////flags
	bool verbose=true;

	////linear solver
	Field<int,d> grid_to_matrix;
	Array<int> matrix_to_grid;
	SparseMatrixT A;
	VectorX p;				////unknown
	VectorX div_u;			////rhs: the solver solves Ap=div_u
	bool is_A_initialized=false;
	bool update_A=true;

	////multigrid solver
	Field<short, d> mat_id;
	bool is_irregular_domain = false;
	MultiGrid::Params multigrid_params;
	GMGPCG_Solver_CPU<d> gmg_solver_cpu;
	#ifdef USE_CUDA
	GMGPCG_Solver_GPU<real,d> gmg_solver_gpu;
	#endif

	////divergence control
	bool use_vol_control = false;
	bool calc_current_vol = true;			////calculate the current vol within apply_vol_control_to_b or not
	real target_vol = (real)-1;				////always needs to be set externally!
	real current_vol = (real)-1;			////needs to be set externally or calculated within apply_vol_control_to_b when calc_current_vol=true			
	real vol_control_ks = (real)1e2;

	// default narrow band width
	int narrow_band_cell_num = 5;
	int dirac_band_cell_num = 3;
	
public:
	////constructors
	ProjectionTwoPhase(MacGrid<d>* _mac_grid, FaceField<real, d>* _velocity, FaceField<real, d>* _rho_face, LevelSet<d>* _levelset, Field<ushort, d>* _type = nullptr, BoundaryConditionMacGrid<d>* _bc = nullptr);
	~ProjectionTwoPhase();

	virtual void Initialize(MacGrid<d>* _mac_grid, FaceField<real, d>* _velocity, FaceField<real, d>* _rho_face, LevelSet<d>* _levelset, Field<ushort, d>* _type, BoundaryConditionMacGrid<d>* _bc);

	////set attributes
	void Set_Velocity(FaceField<real,d>& _velocity){if(velocity!=nullptr&&own_velocity)delete velocity;velocity=&_velocity;own_velocity=false;}
	void Set_Rho_face(FaceField<real,d>& _rho_face){if(rho_face!=nullptr&&own_rho_face)delete rho_face;rho_face=&_rho_face;own_rho_face=false;}
	void Set_Level_Set(LevelSet<d>& _levelset){if(levelset != nullptr && own_levelset)delete levelset; levelset = &_levelset; own_levelset = false;}
	void Set_Type(Field<ushort,d>& _type){if(type!=nullptr&&own_type)delete type;type=&_type;own_type=false;}
	void Set_BC(BoundaryConditionMacGrid<d>& _bc){if(bc!=nullptr&&own_bc)delete bc;bc=&_bc;own_bc=false;}

	////projection functions
	virtual void Allocate_System();
	virtual void Update_A();
	virtual void Update_b();			////calculate b as div velocity
	void Apply_Jump_Condition_To_b();
	void Apply_Vol_Control_To_b();
	void Apply_Implicit_Surface_Tension(const real dt);
	virtual void Correction();
	void Update_Mat_Id();			////for irregular domain
	virtual void Build();					////call allocate, update_A, and update_b
	void Solve_CPX(void);
	virtual void Solve();
	virtual void Project();					////call both build, solve, and correction

	////read data
	void Pressure(Field<real,d>& pressure) const;				////write values of p into pressure
	void Pressure_Gradient(FaceField<real,d>& grad_p) const;	////write values of p into pressure
	void Divergence(Field<real,d>& div) const;					////write values of velocity into div
    inline real Pressure_Jump(const VectorD& pos) const 
	{ real curvature = (*levelset).Curvature(pos); return current_dt * sigma * curvature; }
	
	inline bool Is_Levelset_Interface_Face_Index(const std::pair<int, int> face_idx) const
	{
		int axis = face_idx.first;
		VectorDi face = mac_grid->face_grids[axis].Node_Coord(face_idx.second);
		VectorD pos = mac_grid->Face_Center(axis, face);
		if ((*bc).Is_Psi_N(axis, face)) return false;
		real phi = (*levelset).Phi(pos);
		return (phi > -(real)narrow_band_cell_num * mac_grid->grid.dx && phi < (real)narrow_band_cell_num* mac_grid->grid.dx);
	}

	inline real Levelset_Dirac(const real phi) const
	{
		if (phi < -(real)dirac_band_cell_num * mac_grid->grid.dx) return 0;
		else if (phi > (real)dirac_band_cell_num * mac_grid->grid.dx) return 0;
		else return 0.5 * (1.0 + cos(pi * phi / (real)dirac_band_cell_num / mac_grid->grid.dx)) / ((real)dirac_band_cell_num * mac_grid->grid.dx);
	}

	////Physical interface functions that defines the problem
	//NOTE: if you want to design a derived class of this, theoretically you only need to implement these 5 functions
	virtual real Off_Diag_Term(const VectorDi& fluid_cell, const int& nbidx)const;
	virtual real Diag_Face_Term(const int& axis, const VectorDi& face)const;
	virtual real Velocity_Offset(const int& axis, const VectorDi& face)const;
	virtual bool Is_Valid_Cell(const VectorDi& cell) const {return mac_grid->grid.Valid_Cell(cell)&&!(*bc).Is_Psi_D(cell);}
	virtual bool Is_Fluid_Cell(const VectorDi& cell) const {return Is_Valid_Cell(cell);}
};

#endif
