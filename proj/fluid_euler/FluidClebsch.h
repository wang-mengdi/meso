//////////////////////////////////////////////////////////////////////////
// Clebsch Fluid: Implementation of Schrödinger's Smoke [Chern et al. 2015]
// Copyright (c) (2022-), Bo Zhu, Mengdi Wang, Zhecheng Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once

#include <array>

#include "Multigrid.h"
#include "ConjugateGradient.h"
#include "ExteriorDerivative.h"
#include "Advection.h"
#include "Simulator.h"
#include "IOFunc.h"
#include "AuxFunc.h"
#include "GridEulerFunc.h"
#include "device_launch_parameters.h"
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>

namespace Meso {
	using ArrayDC = thrust::device_vector<C>;
	
	//__device__ void c2_copy(const Vector2& origin, ArrayDC& target) {
	//	//C* targetPtr = thrust::raw_pointer_cast(&target[0]);
	//	for (int i = 0; i < 2; i++) {
	//		//targetPtr[i] = origin[i];
	//	}
	//}
	
	__device__ __host__ C Vector2C_Dot(const Vector2C& a, const Vector2C& b) {
		return a[0] * b[0] + a[1] * b[1];
	}


	__global__ void W2V_Mapping_Kernel2_Padding0(const Grid<2> grid, real* face_x, real* face_y, const Vector2C* cell, const real h_bar_over_dx)
	{
		Typedef_VectorD(2);
		VectorDi coord = GPUFunc::Thread_Coord<2>(blockIdx, threadIdx);

		if (!grid.Valid(coord)) return;
		const Vector2C& cell_data = cell[grid.Index(coord)];
		const Vector2C cell_data_conj(thrust::conj(cell_data[0]), thrust::conj(cell_data[1]));

		int face_ind_x = Neighbor_Face_Index(grid, coord, 0, 1);
		int face_ind_y = Neighbor_Face_Index(grid, coord, 1, 1);

		if (face_ind_x != -1 && grid.Valid(coord + VectorDi::Unit(0))) {
			const Vector2C& nb_cell_data = cell[grid.Index(coord + VectorDi::Unit(0))];
			//face_x[face_ind_x] = h_bar_over_dx * thrust::arg(cell_data_conj.dot(nb_cell_data));
			//const C re_product = ArrayFunc::Conj_Dot(nb_cell_data, cell_data);
			const C re_product = Vector2C_Dot(nb_cell_data, cell_data_conj);
			face_x[face_ind_x] = h_bar_over_dx * thrust::arg(re_product);
		}
		if (face_ind_y != -1 && grid.Valid(coord + VectorDi::Unit(1))) {
			const Vector2C& nb_cell_data = cell[grid.Index(coord + VectorDi::Unit(1))];
			//const C re_product = ArrayFunc::Conj_Dot(nb_cell_data, cell_data);
			//const C re_product = thrust::inner_product(nb_cell_data.begin(), nb_cell_data.end(), cell_data_conj.begin(), C(0.0f, 0.0f));
			const C re_product = Vector2C_Dot(nb_cell_data, cell_data_conj);
			face_y[face_ind_y] = h_bar_over_dx * thrust::arg(re_product);			//Vector2i face_coord_y = grid.Face_Coord(1, face_ind_y);
			//printf("cell_data_conj:\n(%f,%f), (%f,%f)\nnb_cell_data:\n(%f,%f), (%f,%f)\n\nface_coord_y:[%d,%d]\nface_y[face_ind_y]:\n%f\n",
			//	cell_data_conj[0].real(), cell_data_conj[0].imag(), cell_data_conj[1].real(), cell_data_conj[1].imag(),
			//	nb_cell_data[0].real(), nb_cell_data[0].imag(), nb_cell_data[1].real(), nb_cell_data[1].imag(),
			//	face_coord_y[0], face_coord_y[1], face_y[face_ind_y]);
		}
	}

	// for blockDim = (4, 4, 4)
	__global__ void W2V_Mapping_Kernel3_Padding0(const GridIndexer<3> grid, real* face_x, real* face_y, real* face_z, const Vector2C* cell, const real h_bar_over_dx)
	{
		Typedef_VectorD(3);
		VectorDi coord = GPUFunc::Thread_Coord<3>(blockIdx, threadIdx);

		const int idx = threadIdx.x;
		const int idy = threadIdx.y;
		const int idz = threadIdx.z;

		if (!grid.Valid(coord)) return;
		const Vector2C cell_data = cell[grid.Index(coord)];
		//ArrayDC cell_data_d(2);
		//c2_copy(cell_data, cell_data_d);
		//cell_data_d[0] = cell_data[0]; cell_data_d[1] = cell_data[1];
		const Vector2C cell_data_conj(thrust::conj(cell_data[0]), thrust::conj(cell_data[1]));
		//ArrayDC cell_data_conj_d(2);
		//cell_data_conj_d[0] = cell_data_conj[0]; cell_data_conj_d[1] = cell_data_conj[1];
		//c2_copy(cell_data_conj, cell_data_conj_d);

		int face_ind_x = Neighbor_Face_Index(grid, coord, 0, 1);
		int face_ind_y = Neighbor_Face_Index(grid, coord, 1, 1);
		int face_ind_z = Neighbor_Face_Index(grid, coord, 2, 1);

		if (face_ind_x != -1 && grid.Valid(coord + VectorDi::Unit(0))) {
			const Vector2C nb_cell_data = cell[grid.Index(coord + VectorDi::Unit(0))];
			//ArrayDC nb_cell_data_d(2);
			//nb_cell_data_d[0] = nb_cell_data[0]; nb_cell_data_d[1] = nb_cell_data[1];
			//c2_copy(nb_cell_data, nb_cell_data_d);
			//const C re_product = ArrayFunc::Conj_Dot(nb_cell_data, cell_data);
			//const C re_product = thrust::inner_product(nb_cell_data_d.begin(), nb_cell_data_d.end(), cell_data_conj_d.begin(), C(0.0f, 0.0f));
			const C re_product = Vector2C_Dot(nb_cell_data, cell_data_conj);
			face_x[face_ind_x] = h_bar_over_dx * thrust::arg(re_product);
			//Vector3i face_coord_x = grid.Face_Coord(0, face_ind_x);
			//printf("cell_data_conj: (%f,%f), (%f,%f) nb_cell_data: (%f,%f), (%f,%f)  face_coord_x:[%d,%d,%d] face_x[face_ind_x]:\n%f\n",
			//	cell_data_conj[0].real(), cell_data_conj[0].imag(), cell_data_conj[1].real(), cell_data_conj[1].imag(),
			//	nb_cell_data[0].real(), nb_cell_data[0].imag(), nb_cell_data[1].real(), nb_cell_data[1].imag(),
			//	face_coord_x[0], face_coord_x[1], face_coord_x[2], face_x[face_ind_x]);

		}
		if (face_ind_y != -1 && grid.Valid(coord + VectorDi::Unit(1))) {
			const Vector2C nb_cell_data = cell[grid.Index(coord + VectorDi::Unit(1))];
			//ArrayDC nb_cell_data_d(2);
			//nb_cell_data_d[0] = nb_cell_data[0]; nb_cell_data_d[1] = nb_cell_data[1];
			//c2_copy(nb_cell_data, nb_cell_data_d);
			//const C re_product = ArrayFunc::Conj_Dot(nb_cell_data, cell_data);
			//const C re_product = thrust::inner_product(nb_cell_data_d.begin(), nb_cell_data_d.end(), cell_data_conj_d.begin(), C(0.0f, 0.0f));
			const C re_product = Vector2C_Dot(nb_cell_data, cell_data_conj);
			face_y[face_ind_y] = h_bar_over_dx * thrust::arg(re_product);
			//Vector3i face_coord_y = grid.Face_Coord(1, face_ind_y);
			//printf("cell_data_conj: (%f,%f), (%f,%f) nb_cell_data: (%f,%f), (%f,%f)  face_coord_y:[%d,%d,%d] face_y[face_ind_y]:\n%f\n",
			//	cell_data_conj[0].real(), cell_data_conj[0].imag(), cell_data_conj[1].real(), cell_data_conj[1].imag(),
			//	nb_cell_data[0].real(), nb_cell_data[0].imag(), nb_cell_data[1].real(), nb_cell_data[1].imag(),
			//	face_coord_y[0], face_coord_y[1], face_coord_y[2], face_y[face_ind_y]);
		}
		if (face_ind_z != -1 && grid.Valid(coord + VectorDi::Unit(2))) {
			const Vector2C nb_cell_data = cell[grid.Index(coord + VectorDi::Unit(2))];
			//ArrayDC nb_cell_data_d(2);
			//nb_cell_data_d[0] = nb_cell_data[0]; nb_cell_data_d[1] = nb_cell_data[1];
			//c2_copy(nb_cell_data, nb_cell_data_d);
			//const C re_product = ArrayFunc::Conj_Dot(nb_cell_data, cell_data);
			//const C re_product = thrust::inner_product(nb_cell_data_d.begin(), nb_cell_data_d.end(), cell_data_conj_d.begin(), C(0.0f, 0.0f));
			const C re_product = Vector2C_Dot(nb_cell_data, cell_data_conj);
			face_z[face_ind_z] = h_bar_over_dx * thrust::arg(re_product);
			//Info("coord is {}, index is , value is {}", coord + VectorDi::Unit(2), face_ind_z, face_y[face_ind_z]);

		}
	}


	//// Can be optimized by shared memory
	//__global__ void W2V_Mapping_Kernel2_Padding0_new(const Grid<2> grid, real* face_x, real* face_y, const Vector2C* cell, const Vector2C* cell_grad, const real h_bar_over_dx)
	//{
	//	Typedef_VectorD(2);
	//	VectorDi coord = GPUFunc::Thread_Coord<2>(blockIdx, threadIdx);

	//	if (!grid.Valid(coord)) return;
	//	const Vector2C& cell_data = cell[grid.Index(coord)];
	//	// question: cell_grad is face_field but cell is field. How should I ??
	//	Vector2C& cell_grad_data;
	//	const Vector2C cell_data_conj(thrust::conj(cell_data[0]), thrust::conj(cell_data[1]));

	//	int face_ind_x = Neighbor_Face_Index(grid, coord, 0, 1);
	//	int face_ind_y = Neighbor_Face_Index(grid, coord, 1, 1);

	//	if (face_ind_x != -1 && grid.Valid(coord + VectorDi::Unit(0))) {
	//		const Vector2C& cell_grad_data = cell_grad[grid.Index(coord + VectorDi::Unit(0))];
	//		// question : why is thrust::arg here ???
	//		// face_x[face_ind_x] = h_bar_over_dx * thrust::arg(cell_data_conj.dot(nb_cell_data));
	//		face_x[face_ind_x] = h_bar_over_dx * (cell_data_conj.dot(cell_grad_data)).real();
	//	}
	//	if (face_ind_y != -1 && grid.Valid(coord + VectorDi::Unit(1))) {
	//		const Vector2C& cell_grad_data = cell_grad[grid.Index(coord + VectorDi::Unit(1))];
	//		// face_y[face_ind_y] = h_bar_over_dx * thrust::arg(cell_data_conj.dot(nb_cell_data));
	//		face_y[face_ind_y] = h_bar_over_dx * (cell_data_conj.dot(cell_grad_data)).real();

	//	}
	//}

	template<int d>
	__global__ void Wave_Function_Normalization_Kernel(const Grid<d> grid, Vector2C* wave_function) {
		const int index = grid.Index(GPUFunc::Thread_Coord<d>(blockIdx, threadIdx));
		Vector2C psi = wave_function[index];
		real norm = std::sqrt(thrust::norm(psi[0]) + thrust::norm(psi[1]));	// thrust::norm returns norm of magnitude of a complex number
		wave_function[index] = psi / norm;
	}

	template<int d>
	__global__ void Wave_Function_Correction_Kernel(const Grid<d> grid, Vector2C* wave_function, const real* pressure, const real dx_over_h_har) {
		const int index = grid.Index(GPUFunc::Thread_Coord<d>(blockIdx, threadIdx));
		const real cell_p = pressure[index];
		Vector2C psi = wave_function[index];
		C c(thrust::exp(C(0, -1) * cell_p * dx_over_h_har));
		psi[0] *= c;
		psi[1] *= c;
		wave_function[index] = psi;
	}

	template<int d>
	class FluidClebsch : public Simulator {
		Typedef_VectorD(d);
	public:
		real h_bar;
		FieldDv<Vector2C, d> wave_function;
		FaceFieldDv<real, d> velocity;

		BoundaryConditionDirect<FieldDv<Vector2C, d>> psi_D;
		FaceFieldDv<real, d> temp_velocity;
		FieldDv<real, d> pressure;
		FieldDv<real, d> vel_div;
		MaskedPoissonMapping<real, d> poisson;
		VCycleMultigridIntp<real, d> MG_precond;
		ConjugateGradient<real> MGPCG;
		ConjugateGradient<real> CG;
		void Init(real h_bar_, Field<bool, d>& fixed, Field<Vector2C, d>& initial_wave_function, FaceField<real, d>& vol) {
			h_bar = h_bar_;
			wave_function.Deep_Copy(initial_wave_function);
			Exterior_Derivative_W2V(velocity, wave_function);
			//Info("velocity {}", velocity);
			//Info("Fixed {}", fixed);
			//Info("Vol {}", vol);
			psi_D.Init(fixed, initial_wave_function);
			temp_velocity.Init(velocity.grid);
			pressure.Init(velocity.grid);
			vel_div.Init(velocity.grid);

			//poisson.Init(fixed, vol);
			//MG_precond.Init_Poisson(poisson, 2, 2);
			//MGPCG.Init(&poisson, &MG_precond, false, -1, 1e-6);
		}
		virtual real CFL_Time(const real cfl) {
			real dx = velocity.grid.dx;
			real max_vel = GridEulerFunc::Linf_Norm(velocity);
			return dx * cfl / max_vel;
		}
		virtual void Output(DriverMetaData& metadata) {
			std::string vts_name = fmt::format("vts{:04d}.vts", metadata.current_frame);
			bf::path vtk_path = metadata.base_path / bf::path(vts_name);
			VTKFunc::Write_VTS(velocity, vtk_path.string());
		}
		virtual void Schrodinger_Solve() {

		}
		virtual void Normalize() {
			wave_function.grid.Exec_Kernel(
				&Wave_Function_Normalization_Kernel<d>,
				wave_function.grid,
				wave_function.Data_Ptr()
			);
		}
		virtual void Exterior_Derivative_W2V(FaceFieldDv<real, d>& F, const FieldDv<Vector2C, d>& C) {
			Assert(!C.Empty(), "Exterior_Derivative_W2V C->F error: C is empty");
			const real h_bar_over_dx = h_bar / C.grid.dx;
			F.Init(C.grid, MathFunc::Zero<real>());
			const Vector2C* cell = C.Data_Ptr();
			if constexpr (d == 2) C.grid.Exec_Kernel(&W2V_Mapping_Kernel2_Padding0, C.grid, F.Data_Ptr(0), F.Data_Ptr(1), cell, h_bar_over_dx);
			else if constexpr (d == 3) C.grid.Exec_Kernel(&W2V_Mapping_Kernel3_Padding0, C.grid, F.Data_Ptr(0), F.Data_Ptr(1), F.Data_Ptr(2), cell, h_bar_over_dx);
		}
		
		// original transformation method, not used
		//virtual void Exterior_Derivative_W2V_original(FaceFieldDv<real, d>& F, const FieldDv<Vector2C, d>& C) {
		//Assert(!C.Empty(), "Exterior_Derivative_W2V_new C->F error: C is empty");

		//FaceFieldDv<real, d> C_grad;
		//ExteriorDerivativePadding0 gradient;
		//gradient::Apply(C_grad, C);
		//const real h_bar_over_dx = h_bar / C.grid.dx;
		//F.Init(C.grid, MathFunc::Zero<real>());
		//const Vector2C* cell = C.Data_Ptr();
		//if constexpr (d == 2) C.grid.Exec_Kernel(&W2V_Mapping_Kernel2_Padding0_new, C.grid, F.Data_Ptr(0), F.Data_Ptr(1), cell, cell_grad, h_bar_over_dx);
		////else if constexpr (d == 3) C.grid.Exec_Kernel(&W2V_Mapping_Kernel3_Padding0, C.grid, F.Data_Ptr(0), F.Data_Ptr(1), F.Data_Ptr(2), cell, h_bar_over_dx);
		//}

		virtual void Wave_Function_Correction() {
			wave_function.grid.Exec_Kernel(
				&Wave_Function_Correction_Kernel<d>,
				wave_function.grid,
				wave_function.Data_Ptr(),
				pressure.Data_Ptr(),
				wave_function.grid.dx / h_bar
			);
		}
		virtual void Advance(DriverMetaData& metadata) {
			real dt = metadata.dt;

			//TODO: solve 
			Schrodinger_Solve();

			Normalize();

			Info("Vel before W2V {}", velocity);
			//wave function to velocity
			Exterior_Derivative_W2V(velocity, wave_function);
			Info("Vel after W2V {}", velocity);

			//projection
			//vel_div=div(velocity)
			ExteriorDerivativePadding0::Apply(vel_div, velocity);
			Info("Div vel {}", vel_div);

			auto [iter, res] = MGPCG.Solve(pressure.Data(), vel_div.Data());
			Info("Solve poisson with {} iters and residual {}", iter, res);
			Info("Pressure {}", pressure);

			//velocity+=grad(p)
			ExteriorDerivativePadding0::Apply(temp_velocity, pressure);
			temp_velocity *= poisson.vol;

			velocity += temp_velocity;
			Info("Vel after projection {}", velocity);

			ExteriorDerivativePadding0::Apply(vel_div, velocity);
			Info("After projection max div {}", GridEulerFunc::Linf_Norm(vel_div));
			
			Wave_Function_Correction();
			psi_D.Apply(wave_function);
		}
	};
}