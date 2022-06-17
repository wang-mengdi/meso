//////////////////////////////////////////////////////////////////////////
// Discrete Differential Geometry Operators
// Copyright (c) (2022-), Zangyueyang Xian, Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#include "Field.h"
#include "FaceField.h"

#pragma once

namespace Meso {

	template<int d>
	constexpr int __host__ __device__ Neighbor_Face_Index(const GridIndexer<d> grid, Vector<int, d> cell, int axis, int side) {
		//side==0: -
		//side==1: +
		cell[axis] += side;
		return grid.Valid_Face(axis, cell) ? grid.Face_Index(axis, cell) : -1;
	}


	// for blockDim = (8, 8)
	// iterate through cell
	// face_x(i,j)=cell(i+1,j)-cell(i,j)
	// face_y(i,j)=cell(i,j+1)-cell(i,j)
	// seems that it behaves like cell(i,j)=0 for outside the boundary
	template<class T>
	__global__ void D_CoCell_Kernel2_Padding0(const GridIndexer<2> grid, T* face_x, T* face_y, const T* cell)
	{
		Typedef_VectorD(2);
		VectorDi coord = GPUFunc::Thread_Coord<2>(blockIdx, threadIdx);
		const int idx = threadIdx.x;
		const int idy = threadIdx.y;

		//part 1: load data into shared memory
		// 8x8 cell data load
		// cell coord: bx*8+idx, by*8+idy
		__shared__ T shared_cell_data[64];
		shared_cell_data[idy * 8 + idx] = grid.Valid(coord) ? cell[grid.Index(coord)] : MathFunc::Zero<T>();
		__syncthreads();

		const T cell_data = shared_cell_data[idy * 8 + idx];
		const T neg_cell_data = -cell_data;

		int face_ind_x0 = Neighbor_Face_Index(grid, coord, 0, 0);
		int face_ind_x1 = Neighbor_Face_Index(grid, coord, 0, 1);
		int face_ind_y0 = Neighbor_Face_Index(grid, coord, 1, 0);
		int face_ind_y1 = Neighbor_Face_Index(grid, coord, 1, 1);

		//NOTE: we believe we can somehow skip some -1 check here
		//part 2: fill in-block faces
		if (idx != 7 && face_ind_x1 != -1) face_x[face_ind_x1] = shared_cell_data[idy * 8 + (idx + 1)] - cell_data;//x+ face
		if (idy != 7 && face_ind_y1 != -1) face_y[face_ind_y1] = shared_cell_data[(idy + 1) * 8 + idx] - cell_data;//y+ face

		//part 3: accumulate on-boundary faces
		if (idx == 0 && face_ind_x0 != -1) MathFunc::Atomic_Add(face_x + face_ind_x0, cell_data);//x- face
		if (idx == 7 && face_ind_x1 != -1) MathFunc::Atomic_Add(face_x + face_ind_x1, neg_cell_data);//x+ face
		if (idy == 0 && face_ind_y0 != -1) MathFunc::Atomic_Add(face_y + face_ind_y0, cell_data);//y- face
		if (idy == 7 && face_ind_y1 != -1) MathFunc::Atomic_Add(face_y + face_ind_y1, neg_cell_data);//y+ face
	}

	// for blockDim = (4, 4, 4)
	// iterate through cell
	// accumulate cell values with sign on adjacent faces
	template<class T>
	__global__ void D_CoCell_Kernel3_Padding0(const GridIndexer<3> grid, T* face_x, T* face_y, T* face_z, const T* cell)
	{
		Typedef_VectorD(3);
		VectorDi coord = GPUFunc::Thread_Coord<3>(blockIdx, threadIdx);

		const int idx = threadIdx.x;
		const int idy = threadIdx.y;
		const int idz = threadIdx.z;

		__shared__ T shared_cell_data[64];

		// 8x8 cell data load
		shared_cell_data[(idz * 4 + idy) * 4 + idx] = grid.Valid(coord) ? cell[grid.Index(coord)] : MathFunc::Zero<T>();
		__syncthreads();

		const T cell_data = shared_cell_data[(idz * 4 + idy) * 4 + idx];
		const T neg_cell_data = -cell_data;
		int face_ind_x0 = Neighbor_Face_Index(grid, coord, 0, 0);
		int face_ind_x1 = Neighbor_Face_Index(grid, coord, 0, 1);
		int face_ind_y0 = Neighbor_Face_Index(grid, coord, 1, 0);
		int face_ind_y1 = Neighbor_Face_Index(grid, coord, 1, 1);
		int face_ind_z0 = Neighbor_Face_Index(grid, coord, 2, 0);
		int face_ind_z1 = Neighbor_Face_Index(grid, coord, 2, 1);

		//NOTE: we believe we can somehow skip some -1 check here

		if (idx != 3 && face_ind_x1 != -1) face_x[face_ind_x1] = shared_cell_data[(idz * 4 + idy) * 4 + (idx + 1)] - cell_data;//x+ face
		if (idy != 3 && face_ind_y1 != -1) face_y[face_ind_y1] = shared_cell_data[(idz * 4 + (idy + 1)) * 4 + idx] - cell_data;//y+ face
		if (idz != 3 && face_ind_z1 != -1) face_z[face_ind_z1] = shared_cell_data[((idz + 1) * 4 + idy) * 4 + idx] - cell_data;//z+ face

		// x-axis shared faces
		if (idx == 0 && face_ind_x0 != -1) MathFunc::Atomic_Add(face_x + face_ind_x0, cell_data);
		if (idx == 3 && face_ind_x1 != -1) MathFunc::Atomic_Add(face_x + face_ind_x1, neg_cell_data);
		// y-axis shared faces
		if (idy == 0 && face_ind_y0 != -1) MathFunc::Atomic_Add(face_y + face_ind_y0, cell_data);
		if (idy == 3 && face_ind_y1 != -1) MathFunc::Atomic_Add(face_y + face_ind_y1, neg_cell_data);
		// z-axis shared faces
		if (idz == 0 && face_ind_z0 != -1) MathFunc::Atomic_Add(face_z + face_ind_z0, cell_data);
		if (idz == 3 && face_ind_z1 != -1) MathFunc::Atomic_Add(face_z + face_ind_z1, neg_cell_data);
	}

	// for blockDim = (8, 8)
	// iterate through cell
	// velocity in face_x, face_y
	// calculate divergence to cell
	template<class T>
	__global__ void D_Face_Kernel2_Padding0(const GridIndexer<2> grid, T* cell, const T* face_x, const T* face_y)
	{
		//if (!grid.Valid(GPUFunc::Thread_Coord<2>(blockIdx, threadIdx))) return;

		Typedef_VectorD(2);
		VectorDi coord = GPUFunc::Thread_Coord<2>(blockIdx, threadIdx);

		const int idx = threadIdx.x;
		const int idy = threadIdx.y;

		int face_ind_x0 = Neighbor_Face_Index(grid, coord, 0, 0);
		int face_ind_x1 = Neighbor_Face_Index(grid, coord, 0, 1);
		int face_ind_y0 = Neighbor_Face_Index(grid, coord, 1, 0);
		int face_ind_y1 = Neighbor_Face_Index(grid, coord, 1, 1);

		__shared__ T shared_face_data[72];

		//NOTE: we believe we can somehow skip some -1 check here

		T div = 0;
		{
			// 9x8 face_x data load
			shared_face_data[idx * 8 + idy] = face_ind_x0 != -1 ? face_x[face_ind_x0] : 0;
			if (idx == 7) shared_face_data[64 + idy] = face_ind_x1 != -1 ? face_x[face_ind_x1] : 0;
			__syncthreads();
			
			//[idx+1][idy]-[idx][idy]
			div += shared_face_data[(idx + 1) * 8 + idy] - shared_face_data[idx * 8 + idy];
		}
		__syncthreads();

		{
			// 8x9 face_y data load
			shared_face_data[idy * 8 + idx] = face_ind_y0 != -1 ? face_y[face_ind_y0] : 0;
			if (idy == 7) shared_face_data[64 + idx] = face_ind_y0 != -1 ? face_y[face_ind_y1] : 0;
			__syncthreads();

			//[idy+1][idx]-[idy][idx]
			div += shared_face_data[(idy + 1) * 8 + idx] - shared_face_data[idy * 8 + idx];
		}

		if (grid.Valid(coord)) cell[grid.Index(coord)] = div;
	}

	// for blockDim = (4, 4, 4)
	// iterate through cell
	template<class T>
	__global__ void D_Face_Kernel3_Padding0(const GridIndexer<3> grid, T* cell, const T* face_x, const T* face_y, const T* face_z)
	{
		Typedef_VectorD(3);
		VectorDi coord = GPUFunc::Thread_Coord<3>(blockIdx, threadIdx);
		const int idx = threadIdx.x, idy = threadIdx.y, idz = threadIdx.z;
		int face_ind_x0 = Neighbor_Face_Index(grid, coord, 0, 0);
		int face_ind_x1 = Neighbor_Face_Index(grid, coord, 0, 1);
		int face_ind_y0 = Neighbor_Face_Index(grid, coord, 1, 0);
		int face_ind_y1 = Neighbor_Face_Index(grid, coord, 1, 1);
		int face_ind_z0 = Neighbor_Face_Index(grid, coord, 2, 0);
		int face_ind_z1 = Neighbor_Face_Index(grid, coord, 2, 1);

		__shared__ T shared_face_data[80];

		//NOTE: we believe we can somehow skip some -1 check here

		T div = 0;
		{
			// shared addressing order: x-y-z
			// 5x4x4 face_x data load
			shared_face_data[(idx * 4 + idy) * 4 + idz] = face_ind_x0 != -1 ? face_x[face_ind_x0] : 0;
			if (idx == 3) shared_face_data[(4 * 4 + idy) * 4 + idz] = face_ind_x1 != -1 ? face_x[face_ind_x1] : 0;
			__syncthreads();

			// [idx+1][idy][idz]-[idx][idy][idz]
			div += shared_face_data[((idx + 1) * 4 + idy) * 4 + idz] - shared_face_data[(idx * 4 + idy) * 4 + idz];
		}
		__syncthreads();

		{
			// shared addressing order: y-z-x
			// 4x5x4 face_y data load
			shared_face_data[(idy * 4 + idz) * 4 + idx] = face_ind_y0 != -1 ? face_y[face_ind_y0] : 0;
			if (idy == 3) shared_face_data[(4 * 4 + idz) * 4 + idx] = face_ind_y1 != -1 ? face_y[face_ind_y1] : 0;
			__syncthreads();

			// [idy+1][idz][idx]-[idy][idz][idx]
			div += shared_face_data[((idy + 1) * 4 + idz) * 4 + idx] - shared_face_data[(idy * 4 + idz) * 4 + idx];
		}
		__syncthreads();

		{
			// shared addressing order: z-x-y
			// 4x4x5 face_z data load
			shared_face_data[(idz * 4 + idx) * 4 + idy] = face_ind_z0 != -1 ? face_z[face_ind_z0] : 0;
			if (idz == 3) shared_face_data[(4 * 4 + idx) * 4 + idy] = face_ind_z1 != -1 ? face_z[face_ind_z1] : 0;
			__syncthreads();

			// [idz+1][idx][idy]-[idz][idx][idy]
			div += shared_face_data[((idz + 1) * 4 + idx) * 4 + idy] - shared_face_data[(idz * 4 + idx) * 4 + idy];
		}
		__syncthreads();

		if (grid.Valid(coord)) cell[grid.Index(coord)] = div;
	}

	class ExteriorDerivativePadding0 {
	public:
		//input: C is a 0-form on cell (2d or 3d)
		//output: F is a 1-form on face (3d) or 1-form on edge (2d)
		template<class T, int d>
		static void Apply(FaceFieldDv<T, d>& F, const FieldDv<T, d>& C)
		{
			Assert(!C.Empty(), "Exterior_Derivative C->F error: C is empty");
			F.Init(C.grid, MathFunc::Zero<T>());
			//dim3 blocknum, blocksize;
			//auto [block]C.grid.Get_Kernel_Dims(blocknum, blocksize);
			const T* cell = C.Data_Ptr();
			if constexpr (d == 2) C.grid.Exec_Kernel(&D_CoCell_Kernel2_Padding0<T>, C.grid, F.Data_Ptr(0), F.Data_Ptr(1), cell);
			else if constexpr (d == 3) C.grid.Exec_Kernel(&D_CoCell_Kernel3_Padding0<T>, C.grid, F.Data_Ptr(0), F.Data_Ptr(1), F.Data_Ptr(2), cell);
			//if constexpr (d == 2) D_CoCell_Kernel2_Padding0 << <blocknum, blocksize >> > (C.grid, F.Data_Ptr(0), F.Data_Ptr(1), cell);
			//else if constexpr (d == 3) D_CoCell_Kernel3_Padding0 << <blocknum, blocksize >> > (C.grid, F.Data_Ptr(0), F.Data_Ptr(1), F.Data_Ptr(2), cell);
		}

		//input: F is a 2-form on face(3d) or 1-form on edge (2d)
		//output: C is a 3-form on cell (3d) or 2-form on cell (2d)
		template<class T, int d>
		static void Apply(FieldDv<T, d>& C, const FaceFieldDv<T, d>& F) {
			Assert(!F.Empty(), "Exterior_Derivative F->C error: F is empty");
			C.Init(F.grid);
			//dim3 blocknum, blocksize;
			auto [blocknum, blocksize] = F.grid.Get_Kernel_Dims();
			T* cell = C.Data_Ptr();
			if constexpr (d == 2) D_Face_Kernel2_Padding0 << <blocknum, blocksize >> > (F.grid, cell, F.Data_Ptr(0), F.Data_Ptr(1));
			else if constexpr (d == 3) D_Face_Kernel3_Padding0 << <blocknum, blocksize >> > (F.grid, cell, F.Data_Ptr(0), F.Data_Ptr(1), F.Data_Ptr(2));
		}
	};
}