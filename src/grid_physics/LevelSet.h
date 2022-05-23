//////////////////////////////////////////////////////////////////////////
// Level set
// Copyright (c) (2022-), Zhiqi Li
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "Grid.h"
#include "Field.h"
#include "Common.h"
#include "AuxFunc.h"
#include "Interpolation.h"
#include "GeometryPrimitives.h"
#include<iostream>
//#include "cuda_runtime.h"
#include <typeinfo>
/// Some function is called as a intergral and some function is called only for a element
namespace Meso {
	template<int d, class PointIntp, DataHolder side = HOST> class LevelSet {
	public:
		Typedef_VectorD(d);
		Grid<d> grid;
		Grid<d>* grid_ptr;
		Field<real, d, side> phi;
		real* data;
		using Intp = Interpolation<PointIntp>;
	public:
		~LevelSet() {
			if constexpr (side == DEVICE)
				checkCudaErrors(cudaFree(grid_ptr));
		}
		LevelSet() { 
			if constexpr(side==DEVICE)
				checkCudaErrors(cudaMalloc(&grid_ptr, sizeof(Grid<d>)));
		}
		LevelSet(const Grid<d>& _grid):LevelSet() {
			Initialize(_grid);
		}
		__host__ void Initialize(const Grid<d>& _grid) {
			grid = _grid;
			phi = Field<real,d,side>(_grid, std::numeric_limits<real>::max());
			data = thrust::raw_pointer_cast(&((*(phi.data))[0]));
			if constexpr (side == DEVICE)
				checkCudaErrors(cudaMemcpy(grid_ptr, &grid, sizeof(Grid<d>), cudaMemcpyHostToDevice));
			else
				grid_ptr = &grid;
			checkCudaErrors(cudaGetLastError());
			//data = thrust::raw_pointer_cast(phi.data->data());
			auto test = []__device__ __host__(VectorDi coord) {
				return (real)0.0;
			};
			
		}
		template<class PointIntp1, DataHolder side1>
		__host__ void Initialize(LevelSet<d, PointIntp1,side1>& leveset) {
			//Initialize(leveset.grid);
			grid = leveset.grid;
			phi = Field<real, d, side>(grid, std::numeric_limits<real>::max());
			phi.Deep_Copy(leveset.phi);
			data = thrust::raw_pointer_cast(&((*(phi.data))[0]));
			if constexpr (side == DEVICE)
				checkCudaErrors(cudaMemcpy(grid_ptr, &grid, sizeof(Grid<d>), cudaMemcpyHostToDevice));
			else
				grid_ptr = &grid;
			checkCudaErrors(cudaGetLastError());
		}
		
		/// Note!!!!!
		/// This function could only be called for HOST levelset!
		/// This is because ImplicitGeometry have virtual function,
		/// If a object of a class with virtual member function will not crush, 
		/// That object must be malloc in the device!
		__host__ void Set_By_Geom(ImplicitGeometry<d>& geom) {
			/// A strange error: For this host platform, an extended lambda cannot be defined inside the 'if'
			/// or 'else' block of a constexpr if statement
			Assert(side == HOST, "This function could not be called on device!");
			Grid<d>* grid_local_ptr = grid_ptr;
			ImplicitGeometry<d>* geom_ptr = &geom;
			auto calc_f = [geom_ptr, grid_local_ptr]__host__ __device__(const VectorDi& cell) {
				return geom_ptr->Phi(grid_local_ptr->Position(cell));
			};
			phi.Calc_Nodes(calc_f);
			data = thrust::raw_pointer_cast(&((*(phi.data))[0]));
			checkCudaErrors(cudaGetLastError());
		}
		///Note!!!!!!!!!
		/// Following functions, cannot be called from HOST, if side==DEVICE.
		/// This is because, See the Intp:Value function,
		/// Although Intp:Value is host as well, the argument phi is device,
		/// So the statement "const T* data_ptr = F.Data_Ptr();" in Intp:Value will get a pointer on device
		/// while   Intp:Value is a host function, which is a conflict!
		__host__ real Phi(const VectorD& pos) const {
			return Phi(grid_ptr, data, pos);
		}
		__host__ __device__ real Phi(Grid<d>* grid_ptr,real* data,const VectorD& pos) const {
			return Intp::Value(*grid_ptr, data, pos);
		}
		__host__ void All_Phi(Field<real, d, side>& t) const {
			///Note here:!!!!
			/// We cannot capture a data member, which is a pointer to object on device, by this pointer
			/// Assigning it to local pointer and then captured by lambda is required!
			Grid<d>* grid_local_ptr = grid_ptr;
			real* data_local_ptr = data;
			auto calc_f = [grid_local_ptr, data_local_ptr, this]__device__ __host__(VectorDi cell)->real {
				const VectorD& pos = grid_local_ptr->Position(cell);
				//real x = Intp::Value(*grid_local_ptr, data_local_ptr, pos);
				real x = Phi(grid_local_ptr, data_local_ptr, pos);
				return x;
			};
			t.Calc_Nodes(calc_f);
		}


		__host__ VectorD Normal(const VectorD& pos) const {
			////TOFIX: fix finite difference on the boundary
			return Normal(grid_ptr, data, pos);
			/*VectorD normal;
			for (int i = 0; i < d; i++)
				normal[i] = (Phi(pos + VectorD::Unit(i) * grid_ptr->dx) - Phi(pos - VectorD::Unit(i) * grid_ptr->dx)) / ((real)2 * grid_ptr->dx);
			return normal.normalized();*/
		}
		__host__ __device__ VectorD Normal(Grid<d>* grid_ptr, real* data,const VectorD& pos) const {
			////TOFIX: fix finite difference on the boundary
			VectorD normal;
			for (int i = 0; i < d; i++)
				normal[i] = (Phi(grid_ptr,data,pos + VectorD::Unit(i) * grid_ptr->dx) - Phi(grid_ptr, data, pos - VectorD::Unit(i) * grid_ptr->dx)) / ((real)2 * grid_ptr->dx);
			return normal.normalized();
		}
		__host__ void All_Normal(Field<VectorD, d, side>& t) const {
			Grid<d>* grid_local_ptr = grid_ptr;
			real* data_local_ptr = data;
			auto calc_f = [grid_local_ptr, data_local_ptr, this]__device__ __host__(VectorDi cell)->VectorD {
				const VectorD& pos = grid_local_ptr->Position(cell);
				//real x = Intp::Value(*grid_local_ptr, data_local_ptr, pos);
				return Normal(grid_local_ptr, data_local_ptr, pos);
			};
			t.Calc_Nodes(calc_f);
		}
		__host__ void All_Normal(FaceField<real, d, side>& normals) const {
			Grid<d>* grid_local_ptr = grid_ptr;
			real* data_local_ptr = data;
			normals.Calc_Faces(
				[grid_local_ptr, data_local_ptr, this] __device__ __host__(int axis, VectorDi face)->real {
					const VectorD& pos = (grid_local_ptr->Face_Grid(axis)).Position(face);
					return Normal(grid_local_ptr, data_local_ptr, pos)[axis];
				}
			);
		}

		__host__ VectorD Gradient(const VectorD& pos) const {
			////TOFIX: fix finite difference on the boundary
			//VectorD normal;
			//for (int i = 0; i < d; i++)
			//	normal[i] = (Phi(pos + VectorD::Unit(i) * grid_ptr->dx) - Phi(pos - VectorD::Unit(i) * grid_ptr->dx)) / ((real)2 * grid_ptr->dx);
			//return normal;
			return Gradient(grid_ptr, data, pos);
		}
		__host__ __device__ VectorD Gradient(Grid<d>* grid_ptr, real* data,const VectorD& pos) const {
			////TOFIX: fix finite difference on the boundary
			VectorD normal;
			for (int i = 0; i < d; i++)
				normal[i] = (Phi(grid_ptr, data, pos + VectorD::Unit(i) * grid_ptr->dx) - Phi(grid_ptr, data, pos - VectorD::Unit(i) * grid_ptr->dx)) / ((real)2 * grid_ptr->dx);
			return normal;
		}
		__host__ void All_Gradient(Field<VectorD, d, side>& t) const {
			Grid<d>* grid_local_ptr = grid_ptr;
			real* data_local_ptr = data;
			auto calc_f = [grid_local_ptr, data_local_ptr, this]__device__ __host__(VectorDi cell)->VectorD {
				const VectorD& pos = grid_local_ptr->Position(cell);
				//real x = Intp::Value(*grid_local_ptr, data_local_ptr, pos);
				return Gradient(grid_local_ptr, data_local_ptr, pos);
			};
			t.Calc_Nodes(calc_f);
		}

		__host__ VectorD Closest_Point(const VectorD& pos, real epsilon = (real)0) const {
			VectorD normal = Gradient(pos);
			normal.normalize();
			return pos - normal * (Phi(pos) + epsilon);
		}
		__host__ VectorD Closest_Point_With_Iterations(const VectorD& pos, const int max_iter = 5) const {
			VectorD intf_pos = pos;
			for (int i = 0; i < max_iter; i++) {
				intf_pos = Closest_Point(intf_pos);
				if (Phi(intf_pos) < (real)0)return intf_pos;
			}
			return intf_pos;
		}
		__host__ __device__ VectorD Closest_Point(Grid<d>* grid_ptr, real* data,const VectorD& pos, real epsilon = (real)0) const {
			VectorD normal = Gradient(grid_ptr, data, pos);
			normal.normalize();
			return pos - normal * (Phi(grid_ptr,data,pos) + epsilon);
		}
		__host__ __device__ VectorD Closest_Point_With_Iterations(Grid<d>* grid_ptr, real* data,const VectorD& pos, const int max_iter = 5) const {
			VectorD intf_pos = pos;
			for (int i = 0; i < max_iter; i++) {
				intf_pos = Closest_Point(grid_ptr, data, intf_pos);
				if (Phi(grid_ptr, data,intf_pos) < (real)0)return intf_pos;
			}
			return intf_pos;
		}
		
		__host__ __device__ real Curvature(Grid<d>* grid_ptr, real* data, const VectorD& pos) const {
			real one_over_dx = (real)1 / grid_ptr->dx; real one_over_two_dx = (real).5 * one_over_dx; 
			real curvature = (real)0;
			for (int i = 0; i < d; i++) {
				VectorD normal_left = Normal(grid_ptr,data,pos - VectorD::Unit(i) * grid_ptr->dx);
				VectorD normal_right = Normal(grid_ptr, data, pos + VectorD::Unit(i) * grid_ptr->dx);
				curvature += (normal_right[i] - normal_left[i]) * one_over_two_dx;
			}
			return curvature;
//			return abs(curvature) < one_over_dx ? curvature : (curvature <= (real)0 ? (real)-1 : (real)1) * one_over_dx;
		}
		__host__ real Curvature(const VectorD& pos) {
			return Curvature(grid_ptr, data, pos);
		}
		__host__ void All_Curvature(Field<real, d, side>& t) const {
			Grid<d>* grid_local_ptr = grid_ptr;
			real* data_local_ptr = data;
			auto calc_f = [grid_local_ptr, data_local_ptr, this]__device__ __host__(VectorDi cell)->real {
				const VectorD& pos = grid_local_ptr->Position(cell);
				//real x = Intp::Value(*grid_local_ptr, data_local_ptr, pos);
				return Curvature(grid_local_ptr, data_local_ptr, pos);
			};
			t.Calc_Nodes(calc_f);
		}
		
		__host__ __device__ real Cell_Fraction(const VectorDi& cell) const {		////approximate cell volume using phi value
			return (real).5 - MathFunc::Clamp(Phi(grid_ptr->Position(cell)), -(real).5 * grid_ptr->dx, (real).5 * grid_ptr->dx) / grid_ptr->dx;
		}
		__host__ real Total_Volume() const {
			real total_vol = (real)0;
			real cell_vol = std::pow(grid.dx, d);
			if constexpr(side==HOST){
				for (int c = 0; c < grid.DoF(); c++) {
					const VectorDi& cell = grid.Coord(c);
					total_vol += Cell_Fraction(cell);
				}
			}
			else {
				///TODO
			}
			return total_vol * cell_vol;
		}
		
		////Helper functions
		
		__host__ __device__ real Sign(const real phi) { return phi <= (real)0 ? (real)-1 : (real)1; }
		__host__ __device__  bool Interface(const real phi_1, const real phi_2) { return Sign(phi_1) != Sign(phi_2); }
		__host__ __device__  real Theta(const real phi_1, const real phi_2) { return phi_1 / (phi_1 - phi_2); }
		/// grid_ptr is required, or the function could not be called from GPU
		/// Remember: lambda cannot capture data member by capturing this, if this is on HOST
		__host__ __device__  real CorrectPhi(Grid<d>* grid_ptr, const VectorDi& cell, ushort* on_interface) {
			VectorD correct_phi = VectorD::Ones() * std::numeric_limits<real>::max();
			VectorDi correct_axis = VectorDi::Zero();
			for (int j = 0; j < 2 * d; j++) {
				VectorDi nb = grid_ptr->Nb_C(cell, j);
				if (!grid_ptr->Valid(nb)) continue;
				const int nb_idx = grid_ptr->Index(nb);
				if (on_interface[nb_idx] && Interface(Phi(grid_ptr->Position(cell)), Phi(grid_ptr->Position(nb)))) {
					/// Here assume Phi is lenear in a little distance
					/// Then normalise it to dx, which is has the actual unit of distance
					/// c_phi is always non-negative
					real c_phi = Theta(Phi(grid_ptr->Position(cell)), Phi(grid_ptr->Position(nb))) * grid_ptr->dx;
					//printf("%f %f %f\n", Phi(grid_ptr->Position(cell)), Phi(grid_ptr->Position(nb)),c_phi);
					int axis = grid_ptr->Nb_C_Axis(j);
					correct_axis[axis] = 1;
					correct_phi[axis] = std::min(correct_phi[axis], c_phi);
				}
			}
			if (correct_axis != VectorDi::Zero()) {
				real hmnc_mean = (real)0;
				for (int i = 0; i < d; i++) {
					if (correct_axis[i] == 0)continue;
					hmnc_mean += (real)1 / (correct_phi[i] * correct_phi[i]);
				}
				hmnc_mean = sqrt((real)1 / hmnc_mean);
				return hmnc_mean;
			}
			else {
				Assert(true, "Error: [Levelset] bad preconditioning");
			}
		}
		__host__ __device__  ushort OnInterface(Grid<d>* grid_ptr,real* data, const VectorDi& cell) {
			ushort ans = 0;
			for (int j = 0; j < 2 * d; j++) {
				VectorDi nb = grid_ptr->Nb_C(cell, j);
				ans = (grid_ptr->Valid(nb) && Interface(Phi(grid_ptr,data,grid_ptr->Position(cell)), Phi(grid_ptr, data,grid_ptr->Position(nb)))) ? 1 : ans;
			}
			return ans;
		}
		__host__ __device__  ushort IsFront(Grid<d>* grid_ptr,ushort* on_interface, const VectorDi& cell) {
			ushort is_front = 0;
			const int idx= grid_ptr->Index(cell);
			for (int j = 0; j < 2 * d; j++) {
				VectorDi nb= grid_ptr->Nb_C(cell, j);
				if (!grid_ptr->Valid(nb))continue;
				const int nb_idx= grid_ptr->Index(nb);
				is_front = (on_interface[nb_idx]==1&&(!on_interface[idx]==1)) ? 1 : is_front;
			}
			return is_front;
		}
		
		//suppose you have a fluid levelset, and negative_levelset is some solid levelset immersed in it
		/*
		template<class PointIntp1>
		__host__  void Fix_With_Complementary(const LevelSet<d, PointIntp1, side>& negative_levelset) {
			const int cell_num = grid.Number_Of_Cells();
			grid.Exec_Each(
				[&](const VectorDi& cell) {
					if (Phi(cell) > 0) return;
					for (int j = 0; j < Grid<d>::Number_Of_Nb_C(); j++) {
						VectorDi nb = grid_ptr->Nb_C(cell, j);
						if (!grid_ptr->Valid(nb)) continue;
//						if (negative_levelset.Interface(negative_levelset.phi(cell), negative_levelset.phi(nb))) {
	//						data[grid->Index(cell)] = -negative_levelset.phi(cell);
//							return;
//						}
					}});
		}
		*/
		// TODO: how to support FIx_With_Complementary in CUDA?
		// 
		//////////////////////////////////////////////////////////////////////////
		////Fast marching

		void Fast_Marching(const real band_width=-1) {
			real InitValue = band_width < 0 ? std::numeric_limits<real>::max() : band_width;
			Field<real, d, side> tent(grid.Counts(), InitValue);
			ArrayPtr<ushort, side> done = std::make_shared<Array<ushort, side>>(grid.DoF(), 0);

			using PRI = std::pair<real, int>;
			std::priority_queue<PRI, Array<PRI>, std::greater<PRI> > heap;

			thrust::counting_iterator<int> idxfirst(0);
			thrust::counting_iterator<int> idxlast = idxfirst + grid.DoF();
			Grid<d>* grid_local_ptr = grid_ptr;
			ushort* done_local_ptr = thrust::raw_pointer_cast(&((* done)[0]));
			real* tent_local_ptr = thrust::raw_pointer_cast(&(*(tent.data))[0]);
			real* data_local_ptr = data;
			thrust::transform(
				idxfirst,
				idxlast,
				done->begin(),
				[grid_local_ptr, data_local_ptr, this]__device__ __host__(const int idx)->ushort {
					const VectorDi cell = grid_local_ptr->Coord(idx);
					return this->OnInterface(grid_local_ptr, data_local_ptr,cell);
				}
			);
			thrust::transform(
				idxfirst,
				idxlast,
				tent.data->begin(),
				[grid_local_ptr, done_local_ptr, InitValue, data_local_ptr, this]__device__ __host__(const int idx)->real {
					const VectorDi cell = grid_local_ptr->Coord(idx);
					if (!done_local_ptr[idx]) return InitValue;
					else return this->CorrectPhi(grid_local_ptr, cell, done_local_ptr);
				}
			);
			thrust::transform(
				idxfirst,
				idxlast,
				tent.data->begin(),
				[grid_local_ptr, done_local_ptr, tent_local_ptr, InitValue, this]__device__ __host__(const int idx) {
					VectorDi cell = grid_local_ptr->Coord(idx);
					if (IsFront(grid_local_ptr, done_local_ptr, cell)) {
						real temp = Solve_Eikonal(cell, tent_local_ptr, grid_local_ptr, done_local_ptr);
						done_local_ptr[idx] =2;
						return temp;
					}
					else return tent_local_ptr[idx];
				}
			); 

			/// Then we must copy the data to host, because the following action must be on CPU
			Field<real, d, HOST> tent_host(tent.grid);
			ArrayPtr<ushort, HOST> done_host;
			if constexpr (side == HOST) {
				tent_host.data = tent.data;
				done_host = done;
			}
			else {
				*(tent_host.data) = *(tent.data);
				done_host = std::make_shared<Array<ushort, HOST>>();
				*done_host = *done;
			}
			for (int i = 0; i < grid.DoF(); i++) {
				real x;
				x = (*(tent_host.data))[i];
				//printf("**%f ", x);
			}
#pragma omp parallel for
			for (int idx = 0; idx < grid.DoF(); idx++) {
				const VectorDi cell = grid.Coord(idx);
				if ((*done_host)[idx] == 2) {
					(*done_host)[idx] = 1;
//					printf("(%f,(%f %f), %f),", (*(tent_host.data))[idx], grid.Position(grid.Coord(idx))[0], grid.Position(grid.Coord(idx))[1], grid.Position(grid.Coord(idx)).norm()-1);
//					printf("(%f,%f) ", (*(tent_host.data))[idx], grid.Position(grid.Coord(idx)).norm() - 1);

					heap.push(PRI((*(tent_host.data))[idx], idx));
				}
			}
			//// heap traversing
			while (!heap.empty()) {
				const real top_val = heap.top().first;
				const int cell_idx = heap.top().second;
				const VectorDi cell = grid.Coord(cell_idx);
				heap.pop();
				if (((*(tent_host.data))[cell_idx]-top_val)> std::numeric_limits<real>::min()*2)continue;
				(*done_host)[cell_idx] = 1;
				for (int i = 0; i < 2 * d; i++) {
					VectorDi nb = grid.Nb_C(cell, i);
					if (!grid.Valid(nb))continue;
					const int nb_idx = grid.Index(nb);
					if (!(*done_host)[nb_idx]) {
						real temp = Solve_Eikonal(nb, tent_host, *done_host);
						
						if (temp < (*(tent_host.data))[nb_idx]) {
							(*(tent_host.data))[nb_idx] = temp;
							heap.push(PRI(temp, nb_idx));
						}
					}
				}
			}
			for (int i = 0; i < grid.DoF(); i++) {
				real x;
				x = (*(tent_host.data))[i];
			}
			/// Then transfer the tent
			if constexpr (side == DEVICE) *(tent.data) = *(tent_host.data);
			tent_local_ptr = thrust::raw_pointer_cast(&(*(tent.data))[0]);
			data_local_ptr = data;
			thrust::transform(
				idxfirst,
				idxlast,
				phi.data->begin(),
				[data_local_ptr, tent_local_ptr, this]__device__ __host__(const int idx) {
					return tent_local_ptr[idx] * Sign(data_local_ptr[idx]);
				}
			);
		}

		__host__ real Solve_Eikonal(const VectorDi& cell, Field<real, d, HOST>& tent, Array<ushort,HOST>& done) {

			return Solve_Eikonal(cell, 
				thrust::raw_pointer_cast ( & ((*(tent.data))[0])),
				&grid,
				thrust::raw_pointer_cast(&(done[0])));
		}
		__device__ __host__ real Solve_Eikonal(const VectorDi& cell, real* tent_local_ptr, Grid<d>* grid_ptr, ushort* done_local_ptr) {
			// calculate correct phi from nb interface cells
			
			VectorD correct_phi = VectorD::Ones() * std::numeric_limits<real>::max();
			VectorDi correct_axis = VectorDi::Zero();
			for (int i = 0; i < 2 * d; i++) {
				VectorDi nb = grid_ptr->Nb_C(cell, i);
				if (!grid_ptr->Valid(nb)) continue;
				const int nb_idx = grid_ptr->Index(nb);
				if (done_local_ptr[nb_idx]==1) {
					int axis = grid_ptr->Nb_C_Axis(i);
					correct_axis[axis] = 1;
					correct_phi[axis] = std::min(correct_phi[axis], tent_local_ptr[nb_idx]);
				}
			}
			// update phi on the cell
			
			real new_phi=0;
			int n = correct_axis.sum();
			
			switch (n) {
			case 1: {
				real c_phi;
				for (int i = 0; i < d; i++)
					if (correct_axis[i] != 0) { c_phi = correct_phi[i]; break; }
				new_phi = grid_ptr->dx + c_phi;
			} break;
			case 2: {
				real p[2];
				int j = 0;
				for (int i = 0; i < d; i++)
					if (correct_axis[i] != 0) p[j++] = correct_phi[i];
				Solve_Quadratic(p[0], p[1], grid_ptr->dx, new_phi);
			} break;
			case 3: {
				Solve_Quadratic(correct_phi[0], correct_phi[1], correct_phi[2], grid_ptr->dx, new_phi);
			} break;
			default: {
				printf("Error: [Levelset] bad solving Eikonal");
				//std::exit(1);
			} break;
			}
			return new_phi;
		}
		__device__ __host__ bool Solve_Quadratic(const real p1, const real p2, const real dx, real& rst) {
			if (abs(p1) >= abs(p2) + dx) { rst = p2 + dx; return true; }
			else if (abs(p2) >= abs(p1) + dx) { rst = p1 + dx; return true; }
			else {
				real delta = (real)2 * dx * dx - pow(p1 - p2, 2);
				Assert(delta >= (real)0, "Error: [Levelset] negative delta in Solve_Quadratic_2");
				rst = (real).5 * (p1 + p2 + sqrt(delta)); return true;
			}
		}
		__device__ __host__ bool Solve_Quadratic(const real p1, const real p2, const real p3, const real dx, real& rst) {
			real delta = pow(p1 + p2 + p3, 2) - (real)3 * (p1 * p1 + p2 * p2 + p3 * p3 - dx * dx);
			if (delta < (real)0) {
				int i = 0; real p_max = abs(p1); if (abs(p2) > p_max) { i = 1; p_max = abs(p2); }if (abs(p3) > p_max) { i = 2; p_max = abs(p3); }
				real q1, q2; if (i == 0) { q1 = p2; q2 = p3; }
				else if (i == 1) { q1 = p1; q2 = p3; }
				else { q1 = p1; q2 = p2; }
				return Solve_Quadratic(q1, q2, dx, rst);
			}
			rst = (real)1.0/3.0 * (p1 + p2 + p3 + sqrt(delta)); return true;
		}

		/*
		///Code delete for this version
				__host__ void Update_Normal(Field<VectorD, d>& normals, int idx) const {
			const VectorDi& cell = grid->Coord(idx);
			const VectorD& pos = grid->Center(cell);
			(*( normals.data))[grid->Index(cell)] = Normal(pos);
		}
		// This function is for kernel
		template<class T>
		__device__ void Update_Normal(T* data, int idx,Grid<d>* g) const {
			const VectorDi& cell = grid_ptr->Coord(idx);
			const VectorD& pos = grid_ptr->Center(cell);
			data[(*g).Index(cell)]= Normal(pos);
		}
		*/
	};
}