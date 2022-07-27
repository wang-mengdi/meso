#pragma once
#include <array>
#include "AuxFunc.h"
#include "GridEulerFunc.h"
#include "device_launch_parameters.h"
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>

namespace Meso {
	namespace WaveFunc {

		__host__ __device__ C Vector2C_Dot(const Vector2C& a, const Vector2C& b);

		__host__ __device__ Vector4 Vector2C_To_Vector4(const Vector<C, 2>& v);

		__host__ __device__ Vector<C, 2> Vector4_To_Vector2C(const Vector4& v);

		template<int d>
		constexpr int __host__ __device__ Neighbor_Face_Idx(const GridIndexer<d> grid, Vector<int, d> cell, int axis, int side) {
			//side==0: -
			//side==1: +
			cell[axis] += side;
			return grid.Valid_Face(axis, cell) ? grid.Face_Index(axis, cell) : -1;
		}

		__global__ void W2V_Mapping_Kernel2_Padding0(const Grid<2> grid, real* face_x, real* face_y, const Vector2C* cell, const real h_bar_over_dx);

		// for blockDim = (4, 4, 4)
		__global__ void W2V_Mapping_Kernel3_Padding0(const GridIndexer<3> grid, real* face_x, real* face_y, real* face_z, const Vector2C* cell, const real h_bar_over_dx);

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
		Vector2C Vel_To_Psi_C(const Vector<real, d>& vel, const Vector<real, d>& pos) {
			Vector2C psi; psi[0] = C(1., 0.); psi[1] = C(.1, 0.);
			real norm = sqrt(thrust::norm(psi[0]) + thrust::norm(psi[1]));
			psi[0] /= norm; psi[1] /= norm;
			real phase = vel.dot(pos);
			for (int i = 0; i < 2; i++) { psi[i] *= exp(i * phase); }
			return psi;
		}

		template<int d>
		void Exterior_Derivative_W2V(FaceFieldDv<real, d>& F, const FieldDv<Vector2C, d>& C, real h_bar) {
			Assert(!C.Empty(), "Exterior_Derivative_W2V C->F error: C is empty");
			const real h_bar_over_dx = h_bar / C.grid.dx;
			F.Init(C.grid, MathFunc::Zero<real>());
			const Vector2C* cell = C.Data_Ptr();
			if constexpr (d == 2) C.grid.Exec_Kernel(&WaveFunc::W2V_Mapping_Kernel2_Padding0, C.grid, F.Data_Ptr(0), F.Data_Ptr(1), cell, h_bar_over_dx);
			else if constexpr (d == 3) C.grid.Exec_Kernel(&WaveFunc::W2V_Mapping_Kernel3_Padding0, C.grid, F.Data_Ptr(0), F.Data_Ptr(1), F.Data_Ptr(2), cell, h_bar_over_dx);
		}
	}
}