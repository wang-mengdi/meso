//////////////////////////////////////////////////////////////////////////
// Two phase incompressible fluid solver on an Eulerian grid
// Copyright (c) (2021-), Bo Zhu Shiying Xiong Zhecheng Wang
// This file is part of CompleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "MacGrid.h"
#include "Advection.h"
#include "BoundaryCondition.h"
#include "ProjectionTwoPhase.h"
#include "SPX_Timer.h"
#include "ClebschHelpers.h"
#include "SPX_Interpolation.h"
#include "FluidFunc.h"
#include "LevelSet.h"

template<int d> class FluidEulerTwoPhaseClebsch
{Typedef_VectorDii(d);
public:
	MacGrid<d> mac_grid;
	FaceField<real,d> velocity;
	
	Field<ushort, d> type;		////supporting type: Fluid (both air and liquid), Solid, Source
	BoundaryConditionMacGrid<d> bc = BoundaryConditionMacGrid<d>(mac_grid);
	real current_time = (real)0;

	//level set
	LevelSet<d> levelset;
	int narrow_band_cell_num = 3;
	real narrow_band_width;

	//wave functions
	real h_bar = (real)0.1;
	Field<Vector4, d> psi_L;
	Field<Vector4, d> psi_A;
	std::function<Vector4(const VectorD&)> initial_psi;

	//two phase density
	FaceField<real,d> rho_face;		////rho on faces
	real rho_L = (real)1.;
	real rho_A = (real)0.001;
	real sqrt_rho_L;
	real sqrt_rho_A;

    //blending coefficient
    real beta = 0.05;
	real blending_band_width = (real)0;
	
	//projection
	ProjectionTwoPhase<d> projection;
	
	//solver settings
	bool verbose = true;
	bool use_zero_extrapolation = true; //zero extrapolation for pos out of field during advection
	bool use_body_force = false;
	
	VectorD g=VectorD::Unit(1)*(real)-9.8;
	bool use_velocity_field = false;	////"velocity" or "wave"

	//evaluation
	std::conditional_t<d == 2, Field<real, 2>, Field<VectorD, 3> > vorticity;
	Field<real, d> divegence_field;

    FluidEulerTwoPhaseClebsch():projection(&mac_grid,&velocity,&rho_face,&levelset,&type,&bc){}

	virtual void Initialize(const VectorDi& cell_counts,const real dx,const VectorD& domain_min=VectorD::Zero())
	{
		mac_grid.Initialize(cell_counts,dx,domain_min);
		//verbose=false;
		//projection.verbose=false;
		projection.Jump_Condition = std::bind(&FluidEulerTwoPhaseClebsch::Pressure_Jump_On_Interface, this, std::placeholders::_1);
		sqrt_rho_L = std::sqrt(rho_L);
		sqrt_rho_A = std::sqrt(rho_A);
		Initialize_Fields();
		Initialize_Interface();
	}

	virtual void Initialize_Wave_Func()
	{
		int cell_num=mac_grid.grid.Number_Of_Cells();
		#pragma omp parallel for
		for(int i=0;i<cell_num;i++){
			VectorDi cell=mac_grid.grid.Cell_Coord(i);
			psi_L(cell) = initial_psi(mac_grid.grid.Center(cell));
			psi_A(cell) = initial_psi(mac_grid.grid.Center(cell));}
		Normalize();
	}

	virtual void Initialize_Fields()
	{
		rho_face.Resize(mac_grid.grid.cell_counts,(real)0);
		velocity.Resize(mac_grid.grid.cell_counts,(real)0);
		psi_L.Resize(mac_grid.grid.cell_counts);
		psi_A.Resize(mac_grid.grid.cell_counts);
		type.Resize(mac_grid.grid.cell_counts,(ushort)CellType::Fluid);
		vorticity.Resize(mac_grid.grid.cell_counts);
		divegence_field.Resize(mac_grid.grid.cell_counts);
	}

	virtual void Initialize_Interface()
	{
		levelset.Initialize(mac_grid.grid);
		narrow_band_width = mac_grid.grid.dx*(real)narrow_band_cell_num;
	}

	virtual void Advance(const real dt)
	{
        Timer timer;timer.Reset();
		Advection(dt);
		if(verbose)timer.Elapse_And_Output_And_Reset("Advection");
		if(!use_velocity_field)Blend_Velocity();
		if(verbose)timer.Elapse_And_Output_And_Reset("Blend_Velocity");	
		if(verbose){std::cout << "before projection" << std::endl;Divergence_Power();}
		Enforce_Incompressibility(dt);
		if(verbose){std::cout << "after projection" << std::endl;Divergence_Power();}
		if(verbose)timer.Elapse_And_Output_And_Reset("Projection");		
		Extrapolation();
		if(verbose)timer.Elapse_And_Output_And_Reset("Extrapolation");
		current_time += dt;
	}

	virtual void Advection(const real dt)
	{
        ////advect interface
		Field<real,d> ghost_phi=levelset.phi;
		std::function<bool(const VectorDi&)> valid_cell=[&](const VectorDi& cell)->bool{return Is_Fluid_Cell(cell);};
		Advection::Semi_Lagrangian(dt,velocity,mac_grid,ghost_phi,mac_grid,levelset.phi,false,valid_cell);
		levelset.Fast_Marching(narrow_band_width);
        Update_Rho_Face();
		Update_Cell_Types();

		////advect velocity
		MacGrid<d> ghost_grid=mac_grid;FaceField<real,d> ghost_velocity=velocity;
		Advection::Semi_Lagrangian(dt, ghost_velocity, mac_grid, ghost_velocity, mac_grid, velocity, use_zero_extrapolation);
		
		////advect psi
		if(!use_velocity_field){
			Field<Vector4,d> ghost_psi_L=psi_L;
			Field<Vector4,d> ghost_psi_A=psi_A;
			Interpolation<d> intp(mac_grid.grid);
			int cell_num=mac_grid.grid.Number_Of_Cells();
			#pragma omp parallel for
			for(int i=0;i<cell_num;i++){
				VectorDi cell=mac_grid.grid.Cell_Coord(i);
				if(!Is_Fluid_Cell(cell))continue;	
				VectorD pos=mac_grid.grid.Center(cell);
				VectorD vel=intp.Interpolate_Face_Vectors(ghost_velocity,pos);
				VectorD mid_pos=pos-vel*(real).5*dt;
				vel=intp.Interpolate_Face_Vectors(ghost_velocity,mid_pos);
				VectorD backtraced_pos=pos-vel*dt;
				VectorDi backtraced_cell=mac_grid.grid.Cell_Coord(backtraced_pos);
				real phi = levelset.phi(cell);
				if(phi<0){
					Vector4 advected_psi_L;
					advected_psi_L=intp.Interpolate_Centers(ghost_psi_L,backtraced_pos);
					Vector<C,2> psi_c_L = V2C(advected_psi_L);
					C c_L = std::exp(1i*0.5*(vel.dot(vel))*dt/h_bar);
					for (int i = 0; i < 2; i++) {psi_c_L[i]*=c_L;} 
					psi_L(cell)=C2V(psi_c_L).normalized()*sqrt_rho_L;}
				else{
					Vector4 advected_psi_A;
					advected_psi_A=intp.Interpolate_Centers(ghost_psi_A,backtraced_pos);
					Vector<C,2> psi_c_A = V2C(advected_psi_A);
					C c_A = std::exp(1i*0.5*(vel.dot(vel))*dt/h_bar);
					for (int i = 0; i < 2; i++) {psi_c_A[i]*=c_A;} 
					psi_A(cell)=C2V(psi_c_A).normalized()*sqrt_rho_A;}}}
	}

	virtual void Extrapolation()
	{
		if (use_velocity_field) return;
		Interpolation<d> intp_vel(mac_grid);
		Interpolation<d> intp_psi(mac_grid.grid);
		int cell_num = mac_grid.grid.Number_Of_Cells();
		#pragma omp parallel for
		for (int i = 0; i < cell_num; i++) {
			VectorDi cell = mac_grid.grid.Cell_Coord(i);
			VectorD pos = mac_grid.grid.Center(cell);
			real phi=levelset.Phi(pos);
			real epsilon = mac_grid.grid.dx;
			VectorD normal=levelset.Gradient(pos);
			normal.normalize();
			if (phi >= (real)0) {
				if(phi < narrow_band_width){
					VectorD interface_pos =  pos-normal*(phi+epsilon);
					VectorD interface_vel=intp_vel.Interpolate_Face_Vectors(velocity,interface_pos);
					Vector4 interface_psi = intp_psi.Interpolate_Centers(psi_L, interface_pos);
					psi_L(cell) = C2V(Vel_To_Psi_C<d>(interface_vel/h_bar,pos-interface_pos,interface_psi)).normalized()*sqrt_rho_L;}
				else{psi_L(cell) = C2V(Vel_To_Psi_C<d>(VectorD::Zero(), pos)).normalized()*sqrt_rho_L;}}
			else if (phi < (real)0) {
				if(phi > -narrow_band_width){
					VectorD interface_pos =  pos-normal*(phi-epsilon);
					VectorD interface_vel=intp_vel.Interpolate_Face_Vectors(velocity,interface_pos);
					Vector4 interface_psi = intp_psi.Interpolate_Centers(psi_A, interface_pos);
					psi_A(cell) = C2V(Vel_To_Psi_C<d>(interface_vel/h_bar,pos-interface_pos,interface_psi)).normalized()*sqrt_rho_A;}
				else{psi_A(cell) = C2V(Vel_To_Psi_C<d>(VectorD::Zero(), pos)).normalized()*sqrt_rho_A;}}}	
		Enforce_Boundary_Conditions();
	}

	virtual void Update_Cell_Types()
	{
		int cell_num=mac_grid.grid.Number_Of_Cells();
		#pragma omp parallel for
		for(int i=0;i<cell_num;i++){
			const VectorDi& cell=mac_grid.grid.Cell_Coord(i);
			type(cell)=(levelset.phi(cell)<=0)?(ushort)CellType::Liquid:(ushort)CellType::Air;
			if(bc.Is_Psi_D(cell))type(cell)=bc.Psi_D_Type(cell);
		}
	}

	virtual void Update_Rho_Face()
	{
		real band_beta = blending_band_width*blending_band_width;
		for (int axis = 0; axis < d; axis++) {
			int face_num = mac_grid.Number_Of_Faces(axis);
			#pragma omp parallel for
			for (int j = 0; j < face_num; j++) {
				VectorDi face = mac_grid.Face_Coord(axis, j);
				VectorDi cell_lr[2];
				real phi_lr[2];
				real tp[2];
				for (int i = 0; i < 2; i++) {
					cell_lr[i] = MacGrid<d>::Face_Incident_Cell(axis, face, i);
					if (!mac_grid.grid.Valid_Cell(cell_lr[i])) { phi_lr[i] = (real)0.; }
					else { phi_lr[i] = levelset.phi(cell_lr[i]); }
					tp[i] = std::max(phi_lr[i], (real)0.);
				}
				real tp_theta = (tp[0] + tp[1]) / (abs(phi_lr[0]) + abs(phi_lr[1]));
				rho_face(axis, face) = rho_A * tp_theta + rho_L * ((real)1. - tp_theta);
				VectorD pos = mac_grid.Face_Center(axis, face); real phi = levelset.Phi(pos);
			}
		}
	}

	void Enforce_Incompressibility(const real dt)
    {
		Enforce_Boundary_Conditions();
		////project velocity
		projection.current_dt=dt;
		projection.Project();
		if (!use_velocity_field) {
			int cell_num = mac_grid.grid.Number_Of_Cells();
			real one_over_dx = (real)1 / mac_grid.grid.dx;
#pragma omp parallel for
			for (int i = 0; i < cell_num; i++) {
				VectorDi cell = mac_grid.grid.Cell_Coord(i);
				if (Is_Fluid_Cell(cell)) {
					real phi0 = levelset.phi(cell);
					//real cell_p = projection.p[projection.grid_to_matrix(cell)];
					real cell_p = projection.meso_pressure_host(cell);
					if (phi0 < 0) {
						Vector<C, 2> psi_c = V2C(psi_L(cell));
						cell_p /= rho_L;
						C c = std::exp(-1i * cell_p * mac_grid.grid.dx / h_bar);
						for (int i = 0; i < 2; i++) { psi_c[i] *= c; }
						psi_L(cell) = C2V(psi_c);
					}
					else {
						Vector<C, 2> psi_c = V2C(psi_A(cell));
						cell_p /= rho_A;
						C c = std::exp(-1i * cell_p * mac_grid.grid.dx / h_bar);
						for (int i = 0; i < 2; i++) { psi_c[i] *= c; }
						psi_A(cell) = C2V(psi_c);
					}
				}
			}
		}
	}

	void Blend_Velocity()
	{
		for (int axis = 0; axis < d; axis++) {
			int face_num = mac_grid.face_grids[axis].Number_Of_Nodes();
			#pragma omp parallel for
			for(int i = 0; i < face_num; i++) {
				VectorDi face = mac_grid.face_grids[axis].Node_Coord(i);
				if(bc.Is_Psi_N(axis,face)) { continue; } 
				const VectorDi cell_0 = mac_grid.Face_Incident_Cell(axis,face,0);
				const VectorDi cell_1 = mac_grid.Face_Incident_Cell(axis,face,1);
                if (!mac_grid.grid.Valid_Cell(cell_0) || !mac_grid.grid.Valid_Cell(cell_1)) {continue;}
				real phi0 = levelset.phi(cell_0);
				real phi1 = levelset.phi(cell_1);
				if (phi0<-blending_band_width && phi1<-blending_band_width) {
					Vector<C, 2> psi_0 = V2C(psi_L(cell_0));
					Vector<C, 2> psi_1 = V2C(psi_L(cell_1));
					C q = psi_0.dot(psi_1);
					velocity(axis, face) = beta*velocity(axis, face) + (1-beta)*std::arg(q)*h_bar/(mac_grid.grid.dx);}
				else if(phi0>blending_band_width && phi1>blending_band_width) {
					Vector<C, 2> psi_0 = V2C(psi_A(cell_0));
					Vector<C, 2> psi_1 = V2C(psi_A(cell_1));
					C q = psi_0.dot(psi_1);
					velocity(axis, face) = beta*velocity(axis, face) + (1-beta)*std::arg(q)*h_bar/(mac_grid.grid.dx);}}}
		Enforce_Boundary_Conditions();
	}
	
	virtual void Enforce_Boundary_Conditions()
	{
		for(auto p:bc.psi_N_values){int axis=p.first[0];int face_index=p.first[1];real value=p.second;
			velocity.face_fields[axis].array[face_index]=value;}
	}
	
	real Pressure_Jump_On_Interface(const VectorD& pos)
	{
		real jump = projection.use_implicit_surface_tension? (real) 0 : projection.sigma * levelset.Curvature(pos);
		//real jump = projection.sigma * levelset.Curvature(pos);
		if (use_body_force) jump -= (rho_L-rho_A) * (g.dot(pos));
		return jump * projection.current_dt;
	}
 
	//////////////////////////////////////////////////////////////////////////
	////helper functions
	real Max_Abs(const FaceField<real,d>& field_q) const
	{
		real max_abs=(real)0;
		for(int i=0;i<d;i++){
			int n=(int)field_q.face_fields[i].array.size();
			for(int j=0;j<n;j++){
				real abs_v=abs(field_q.face_fields[i].array[j]);
				if(abs_v>max_abs)max_abs=abs_v;}}
		return max_abs;
	}

	real Kinetic_Energy(const FaceField<real, d>& field_q) const 
	{
		real kinetic_energy = 0;
		int cell_num = field_q.mac_grid.grid.cell_counts.prod();
		#pragma omp parallel for reduction(+:kinetic_energy)
		for (int i = 0; i < cell_num; i++) {
			const VectorDi cell = mac_grid.grid.Cell_Coord(i);
			VectorD cell_v = VectorD::Zero();
			for (int axis = 0; axis < d; axis++) {
				VectorDi left_face = mac_grid.Cell_Left_Face(axis, cell);
				VectorDi right_face = mac_grid.Cell_Right_Face(axis, cell);
				cell_v[axis] = (real)0.5 * (field_q(axis, left_face) + field_q(axis, right_face));
			}
			kinetic_energy += cell_v.dot(cell_v);
		}
		return (real)0.5 * kinetic_energy * pow(mac_grid.grid.dx, d);
	}

	real Enstrophy()
	{
		FluidFunc::Curl_On_Cell(mac_grid, velocity, vorticity);
		real enstrophy = 0;
		int cell_num = mac_grid.grid.cell_counts.prod();
		#pragma omp parallel for reduction(+:enstrophy)
		for (int i = 0; i < cell_num; i++) {
			const VectorDi cell = mac_grid.grid.Cell_Coord(i);
			if constexpr (d == 2) { enstrophy += pow(vorticity(cell), 2); }
			else { enstrophy += vorticity(cell).dot(vorticity(cell)); }}

		return (real)0.5 * enstrophy * pow(mac_grid.grid.dx, d);
	}

	real Helicity()
	{
		if constexpr (d == 2) return 0.;
		FluidFunc::Curl_On_Cell(mac_grid, velocity, vorticity);
		Field<Vector<real, d>, d> cell_velocity;
		Face_To_Cell_Conversion(velocity, cell_velocity);
		real helicity = 0;
		int cell_num = mac_grid.grid.cell_counts.prod();
		#pragma omp parallel for reduction(+:helicity)
		for (int i = 0; i < cell_num; i++) {
			const VectorDi cell = mac_grid.grid.Cell_Coord(i);
			if constexpr (d == 3) helicity += cell_velocity(cell).dot(vorticity(cell)); }
		return helicity * pow(mac_grid.grid.dx, d);
	}

	void Divergence_Power() {
		real total_divergence = 0.;
		int cell_num = mac_grid.grid.cell_counts.prod();
		#pragma omp parallel for reduction(+:total_divergence)
		for (int i = 0; i < cell_num; i++) {
			VectorDi cell = mac_grid.grid.Cell_Coord(i);
			real divergence = 0.;
			for (int axis = 0; axis < d; axis++) {
			    divergence += velocity(axis, cell + VectorDi::Unit(axis)) - velocity(axis, cell);}
			divergence /= mac_grid.grid.dx;
			double pow_div = pow(divergence, 8);
			total_divergence += pow_div;
			divegence_field(cell) = pow_div;}

		std::cout << "Mean pow 8 of divergence is: " << total_divergence / (real)mac_grid.grid.cell_counts.prod() << std::endl;
	}

	void Normalize()
	{
		int cell_num = mac_grid.grid.Number_Of_Cells();
		#pragma omp parallel for
		for (int i = 0; i < cell_num; i++) {
			VectorDi cell = mac_grid.grid.Cell_Coord(i);
			psi_L(cell) = psi_L(cell).normalized()*sqrt_rho_L;
			psi_A(cell) = psi_A(cell).normalized()*sqrt_rho_A;}
	}
	
    //// cell helpers
	inline bool Is_Fluid_Cell(const VectorDi& cell) const 
	{return mac_grid.grid.Valid_Cell(cell)&&!(type(cell)==(ushort)CellType::Solid||type(cell)==(ushort)CellType::Source);}
	inline bool Is_Solid_Cell(const VectorDi& cell) const {return type(cell)==(ushort)CellType::Solid;}
	inline bool Is_Source_Cell(const VectorDi& cell) const {return type(cell)==(ushort)CellType::Source;}
	inline bool Is_Fluid_Face(const int axis,const VectorDi& face) const
	{{VectorDi cell=MacGrid<d>::Face_Incident_Cell(axis,face,0);if(Is_Fluid_Cell(cell))return true;}
	{VectorDi cell=MacGrid<d>::Face_Incident_Cell(axis,face,1);if(Is_Fluid_Cell(cell))return true;}return false;}
};