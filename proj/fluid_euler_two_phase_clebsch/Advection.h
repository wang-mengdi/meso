//////////////////////////////////////////////////////////////////////////
// Advection schemes (semi-Lagrangian)
// Copyright (c) (2018-), Bo Zhu, Mengdi Wang
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __Advection_h__
#define __Advection_h__
#include "SPX_Field.h"
#include "SPX_FaceField.h"
#include "SPX_Interpolation.h"
#include "SPX_AuxFunc.h"
#include <limits.h>
#include "BoundaryCondition.h"

namespace Advection
{
	template<int d> 
	Vector<real, d> Runge_Kutta2(const real dt, const Interpolation<d>& intp, const FaceField<real, d>& velocity, const Vector<real, d> pos0) {
		Typedef_VectorDii(d);
		VectorD vec0 = intp.Interpolate_Face_Vectors(velocity, pos0);
		VectorD pos1 = pos0 + vec0 * 0.5 * dt;
		VectorD vec1 = intp.Interpolate_Face_Vectors(velocity, pos1);
		return pos0 + vec1 * dt;
	}

	template<int d>
	Vector<real, d> Semi_Lagrangian_Pos(const real dt, std::function<Vector<real, d>(Vector<real, d>)> velocity, const Vector<real, d> pos0) {
		Vector<real, d> vel0 = velocity(pos0);
		Vector<real, d> pos1 = pos0 - vel0 * 0.5 * dt;
		Vector<real, d> vel1 = velocity(pos1);
		return pos0 - vel1 * dt;
	}

////advect a velocity field
template<int d> void Semi_Lagrangian(const real dt,
	const MacGrid<d>& ghost_mac_grid,const FaceField<real,d>& ghost_velocity,
	const MacGrid<d>& mac_grid,FaceField<real,d>& velocity,
	std::function<bool(const int,const Vector<int,d>&)> valid=nullptr)
{	Typedef_VectorDii(d);
	Interpolation<d> ghost_intp(ghost_mac_grid);
	auto vfunc = [&](const VectorD& pos) {return ghost_intp.Interpolate_Face_Vectors(ghost_velocity, pos);};
	//std::bind(&Interpolation<d>::Interpolate_Face_Vectors, &ghost_intp, ghost_velocity, std::placeholders::_1);
	for (int axis = 0;axis < d;axis++) {
		int face_num = mac_grid.Number_Of_Faces(axis);
#pragma omp parallel for
		for (int i = 0;i < face_num;i++) {
			VectorDi face = mac_grid.Face_Coord(axis, i);
			if (valid && !valid(axis, face))continue;	////skip the invalid faces

			Vector<real, d> pos = mac_grid.face_grids[axis].Node(face);			
			Vector<real,d> backtraced_pos = Semi_Lagrangian_Pos<d>(dt,vfunc,pos);
			//VectorD vel = ghost_intp.Interpolate_Face_Vectors(ghost_velocity, pos);
			//VectorD mid_pos = pos - vel * (real).5 * dt;
			//vel = ghost_intp.Interpolate_Face_Vectors(ghost_velocity, mid_pos);
			//VectorD backtraced_pos = pos - vel * dt;
			real advected_v = ghost_intp.Interpolate_Faces(ghost_velocity, backtraced_pos, axis);
			velocity(axis, face) = advected_v;
		}
	}
}

////advect a vector field with referred vector values
////particle position is calculated via ghost_velocity
////old vector value is calculated via ghost_vector
template<int d> void Semi_Lagrangian(const real dt,const FaceField<real,d>& velocity,	/*velocity is on ghost_mac_grid*/
	const MacGrid<d>& ghost_mac_grid,const FaceField<real,d>& ghost_vector,
	const MacGrid<d>& mac_grid,FaceField<real,d>& vector,
	bool use_zero_extrapolation=false,std::function<bool(const int,const Vector<int,d>&)> valid=nullptr)
{	Typedef_VectorDii(d);
	Interpolation<d> intp(ghost_mac_grid);
	std::function<real(const VectorD&, const int)> Ghost_Component = [=](const VectorD& pos, const int axis)->real
	{return (real)0; };

	std::function<VectorD(const VectorD&)> vfunc = nullptr;
	if (use_zero_extrapolation) { vfunc = [&](const VectorD& pos) {return intp.Interpolate_Face_Vectors_With_Ghost(velocity, pos, Ghost_Component);}; }
	else { vfunc = [&](const VectorD& pos) {return intp.Interpolate_Face_Vectors(velocity, pos);}; }
	for(int axis=0;axis<d;axis++){int face_num=mac_grid.Number_Of_Faces(axis);
		#pragma omp parallel for
		for(int i=0;i<face_num;i++){VectorDi face=mac_grid.Face_Coord(axis,i);
			if(valid&&!valid(axis,face))continue;	////skip the invalid faces

			VectorD pos=mac_grid.face_grids[axis].Node(face);
			VectorD backtraced_pos = Semi_Lagrangian_Pos<d>(dt, vfunc, pos);
			//VectorD vel;
			//
			//if (use_zero_extrapolation) { vel=intp.Interpolate_Face_Vectors_With_Ghost(velocity, pos, Ghost_Component); }
			//else{ vel=intp.Interpolate_Face_Vectors(velocity, pos); }
			//VectorD mid_pos=pos-vel*(real).5*dt;
			//if (use_zero_extrapolation) { vel=intp.Interpolate_Face_Vectors_With_Ghost(velocity, mid_pos, Ghost_Component); }
			//else { vel=intp.Interpolate_Face_Vectors(velocity, mid_pos); }
			//VectorD backtraced_pos=pos-vel*dt;
			if (use_zero_extrapolation) { vector(axis, face) = intp.Interpolate_Faces_With_Ghost_Linear(ghost_vector, backtraced_pos, Ghost_Component, axis); }
			else{ vector(axis, face) = intp.Interpolate_Faces(ghost_vector, backtraced_pos, axis); }
		}
	}
}

////advect a cell-based field
template<class T,int d> void Semi_Lagrangian(const real dt,const FaceField<real,d>& velocity,	/*velocity is on ghost_mac_grid*/
	const MacGrid<d>& ghost_mac_grid,const Field<T,d>& ghost_density,
	const MacGrid<d>& mac_grid,Field<T,d>& density,
	bool use_zero_extrapolation=false,std::function<bool(const Vector<int,d>&)> valid=nullptr)
{	Typedef_VectorDii(d);
	Interpolation<d> intp(ghost_mac_grid);
	auto vfunc = [&](const VectorD& pos) {return intp.Interpolate_Face_Vectors(velocity, pos);};

	int cell_num=mac_grid.grid.Number_Of_Cells();
	std::function<real(const VectorD&, const int)> Ghost_Component = [=](const VectorD& pos, const int axis)->real
	{return (real)0; };
	#pragma omp parallel for
	for(int i=0;i<cell_num;i++){VectorDi cell=mac_grid.grid.Cell_Coord(i);
		if(valid&&!valid(cell))continue;		////skip invalid cells

		VectorD pos=mac_grid.grid.Center(cell);
		VectorD backtraced_pos = Semi_Lagrangian_Pos<d>(dt, vfunc, pos);
		//VectorD vel=intp.Interpolate_Face_Vectors(velocity,pos);
		//VectorD mid_pos=pos-vel*(real).5*dt;
		//vel=intp.Interpolate_Face_Vectors(velocity,mid_pos);
		//VectorD backtraced_pos=pos-vel*dt;
			
		T advected_density = Zero<T>();
		if (use_zero_extrapolation) { advected_density = intp.Interpolate_Centers_With_Ghost_Linear(ghost_density, backtraced_pos, Ghost_Component); }
		else { advected_density = intp.Interpolate_Centers(ghost_density, backtraced_pos); }
		density(cell)=advected_density;}
}

////this function only works for scalar advection
template<int d> void Semi_Lagrangian_Catmull_Rom(const real dt,const MacGrid<d>& ghost_grid,const FaceField<real,d>& ghost_velocity,const Field<real,d>& ghost_density,
	const MacGrid<d>& mac_grid,Field<real,d>& density)
{	Typedef_VectorDii(d);
	Interpolation<d> intp(ghost_grid);
	auto vfunc = [&](const VectorD& pos) {return intp.Interpolate_Face_Vectors(ghost_velocity, pos);};
	int cell_num=mac_grid.grid.Number_Of_Cells();
	#pragma omp parallel for
	for(int i=0;i<cell_num;i++){
		VectorDi cell=mac_grid.grid.Cell_Coord(i);
		VectorD pos=mac_grid.grid.Center(cell);
		VectorD backtraced_pos = Semi_Lagrangian_Pos(dt, vfunc, pos);
		//VectorD vel=intp.Interpolate_Face_Vectors(ghost_velocity,pos);
		//VectorD mid_pos=pos-vel*(real).5*dt;
		//vel=intp.Interpolate_Face_Vectors(ghost_velocity,mid_pos);
		//VectorD backtraced_pos=pos-vel*dt;
		real advected_density=intp.Interpolate_Centers_Catmull_Rom(ghost_density,backtraced_pos);
		density(cell)=advected_density;}
}

////advect a vector field
template<int d> void MacCormack(const real dt,const MacGrid<d>& mac_grid,const FaceField<real,d>& velocity,FaceField<real,d>& vector, 
	bool use_zero_extrapolation = false, bool use_clamp=true, BoundaryConditionMacGrid<d>* bc=nullptr)
{	Typedef_VectorDii(d);
	FaceField<real,d> aux_vector=vector;

	Semi_Lagrangian<d>(dt,velocity,mac_grid, velocity,mac_grid,aux_vector,use_zero_extrapolation);
	Semi_Lagrangian<d>(-dt, aux_vector,mac_grid,aux_vector,mac_grid,vector,use_zero_extrapolation);

	if (bc) {
		for (auto p : bc->psi_N_values) {
			int axis = p.first[0]; int face_index = p.first[1]; real value = p.second;
			aux_vector.face_fields[axis].array[face_index] = value;
			vector.face_fields[axis].array[face_index] = value;
		}
	}
	
	for(int axis=0;axis<d;axis++){int face_num=mac_grid.Number_Of_Faces(axis);
		#pragma omp parallel for
		for (int i = 0; i < face_num; i++) {
			VectorDi face = mac_grid.Face_Coord(axis, i);
			if (bc) { // degenerate to semi-lagrangian for cells close to the neumann boundary
				bool neighbor_to_psi_N = false;
				for (int j = 0; j < Grid<d>::Number_Of_Nb_R(); j++) {
					if (bc->Is_Psi_N(axis, mac_grid.face_grids[axis].Nb_R(face, j))) { neighbor_to_psi_N = true; break; }
				}

				if (bc->Is_Psi_N(axis, face)) { continue; }
				else if (neighbor_to_psi_N) { vector(axis, face) = aux_vector(axis, face); }
				else { vector(axis, face) = aux_vector(axis, face) + (real).5 * (velocity(axis, face) - vector(axis, face)); }
			}
			else {
				vector(axis, face) = aux_vector(axis, face) + (real).5 * (velocity(axis, face) - vector(axis, face));
			}
		}
	}
	
	if (use_clamp) {
		Interpolation<d> intp(mac_grid);
		for (int axis = 0; axis < d; axis++) {
			int face_num = mac_grid.Number_Of_Faces(axis);
			#pragma omp parallel for
			for (int i = 0; i < face_num; i++) {
				if (bc&&bc->Is_Psi_N(axis, mac_grid.Face_Coord(axis, i))) { continue; }
				VectorDi face = mac_grid.Face_Coord(axis, i);
				VectorD pos = mac_grid.Face_Center(axis, face);
				VectorD vel = intp.Interpolate_Face_Vectors(velocity, pos);
				VectorD mid_pos = pos - vel * (real).5 * dt;
				vel = intp.Interpolate_Face_Vectors(velocity, mid_pos);
				VectorD backtraced_pos = pos - vel * dt;
				VectorDi backtraced_cell=mac_grid.face_grids[axis].Clamped_Cell_Coord(backtraced_pos); //cell in the face grid
				real min = std::numeric_limits<real>::max();
				real max = -std::numeric_limits<real>::max();
				for (int j = 0; j < Grid<d>::Number_Of_Cell_Incident_Nodes(); j++) {
					VectorDi node = Grid<d>::Cell_Incident_Node(backtraced_cell, j);
					if (vector(axis, node) < min) { min = vector(axis, node); }
					if (vector(axis, node) > max) { max = vector(axis, node); }}
				if (vector(axis, face) < min || vector(axis, face) > max) {//check if overshoot
					vector(axis, face) = intp.Interpolate_Faces(velocity, backtraced_pos, axis);}}}}
}

////advect a scalar field
template<class T, int d> void MacCormack(const real dt, const MacGrid<d>& mac_grid, const FaceField<real, d>& velocity, 
	Field<T, d>& density, bool use_zero_extrapolation = false, bool use_clamp = false)
{
	Typedef_VectorDii(d);
	Field<T, d> ghost_density = density;
	Field<T, d> intermediate_density = density;
	Semi_Lagrangian<T, d>(dt, velocity, mac_grid, ghost_density, mac_grid, intermediate_density,use_zero_extrapolation);
	Semi_Lagrangian<T, d>(-dt, velocity, mac_grid, intermediate_density, mac_grid, density,use_zero_extrapolation);
	int cell_num = mac_grid.grid.Number_Of_Cells();
	#pragma omp parallel for
	for (int i = 0; i < cell_num; i++) {
		VectorDi cell = mac_grid.grid.Cell_Coord(i);
		density(cell) = intermediate_density(cell) + (real).5 * (ghost_density(cell) - density(cell));}

	if (use_clamp) {
		Interpolation<d> intp(mac_grid);
		int cell_num = mac_grid.grid.Number_Of_Cells();
		Grid<d> node_grid=mac_grid.grid.Cell_Grid_To_Node_Grid();
		#pragma omp parallel for
		for (int i = 0; i < cell_num; i++) {
			VectorDi cell = mac_grid.grid.Cell_Coord(i);
			VectorD pos = mac_grid.grid.Center(cell);
			VectorD vel = intp.Interpolate_Face_Vectors(velocity, pos);
			VectorD mid_pos = pos - vel * (real).5 * dt;
			vel = intp.Interpolate_Face_Vectors(velocity, mid_pos);
			VectorD backtraced_pos = pos - vel * dt;
			VectorDi backtraced_node = node_grid.Clamped_Cell_Coord(backtraced_pos);
			real min = std::numeric_limits<real>::max();
			real max = -std::numeric_limits<real>::max();
			for (int j = 0; j < Grid<d>::Number_Of_Cell_Incident_Nodes(); j++) {
				VectorDi node = Grid<d>::Cell_Incident_Node(backtraced_node, j);
				if (density(node) < min) { min = density(node); }
				if (density(node) > max) { max = density(node); }}
			if (density(cell) < min || density(cell) > max) {//check if overshoot
				density(cell) = intp.Interpolate_Centers(ghost_density, backtraced_pos);}}}
}
};

#endif
