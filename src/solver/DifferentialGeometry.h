//////////////////////////////////////////////////////////////////////////
// Discrete Differential Geometry Operators
// Copyright (c) (2022-), Zangyueyang Xian, Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#include "Field.h"
#include "FaceField.h"

namespace Meso {
	// for blockDim = (8, 8)
	// iterate through cell
	// face_x(i,j)=cell(i+1,j)-cell(i,j)
	// face_y(i,j)=cell(i,j+1)-cell(i,j)
	// seems that it behaves like cell(i,j)=0 for outside the boundary
	template<class T>
	__global__ void D_CoCell_Mapping_Kernel2(Grid<2> grid, T* face_x, T* face_y, const T* cell)
	{
		const int nbx = gridDim.x;
		const int nby = gridDim.y;

		const int bx = blockIdx.x;
		const int by = blockIdx.y;

		const int idx = threadIdx.x;
		const int idy = threadIdx.y;

		__shared__ T shared_cell_data[64];

		// 8x8 cell data load
		// cell coord: bx*8+idx, by*8+idy
		// like: shared_cell_data[idx,idy]=cell[cell_coord]
		shared_cell_data[idy * 8 + idx] = cell[grid.Index(Vector2i(bx * 8 + idx, by * 8 + idy))];
		__syncthreads();

		const T cell_data = shared_cell_data[idy * 8 + idx];
		int face_ind;

		// x-axis faces
		face_ind = grid.Face_Index(0, Vector2i(bx * 8 + idx + 1, by * 8 + idy));
		if (idx != 7) face_x[face_ind] = shared_cell_data[idy * 8 + (idx + 1)] - cell_data;

		// y-axis faces
		face_ind = grid.Face_Index(1, Vector2i(bx * 8 + idx, by * 8 + idy + 1));
		if (idy != 7) face_y[face_ind] = shared_cell_data[(idy + 1) * 8 + idx] - cell_data;

		// x-axis shared faces
		face_ind = grid.Face_Index(0, Vector2i(bx * 8 + 0, by * 8 + idy));
		if (idx == 0) atomicAdd(face_x + face_ind, cell_data);

		face_ind = grid.Face_Index(0, Vector2i((bx + 1) * 8 + 0, by * 8 + idy));
		if (idx == 7) atomicAdd(face_x + face_ind, -cell_data);

		// y-axis shared faces
		face_ind = grid.Face_Index(1, Vector2i(bx * 8 + idx, by * 8 + 0));
		if (idy == 0) atomicAdd(face_y + face_ind, cell_data);

		face_ind = grid.Face_Index(1, Vector2i(bx * 8 + idx, (by + 1) * 8 + 0));
		if (idy == 7) atomicAdd(face_y + face_ind, -cell_data);
	}

	// for blockDim = (4, 4, 4)
	// iterate through cell
	template<class T>
	__global__ void D_CoCell_Mapping_Kernel3(Grid<3> grid, T* face_x, T* face_y, T* face_z, const T* cell)
	{
		const int nbx = gridDim.x;
		const int nby = gridDim.y;
		const int nbz = gridDim.z;

		const int bx = blockIdx.x;
		const int by = blockIdx.y;
		const int bz = blockIdx.z;

		const int idx = threadIdx.x;
		const int idy = threadIdx.y;
		const int idz = threadIdx.z;

		__shared__ T shared_cell_data[64];

		// 8x8 cell data load
		shared_cell_data[(idz * 4 + idy) * 4 + idx] = cell[grid.Index(Vector3i(bx * 4 + idx, by * 4 + idy, bz * 4 + idz))];
		__syncthreads();

		const T cell_data = shared_cell_data[(idz * 4 + idy) * 4 + idx];
		int face_ind;

		// x-axis faces
		face_ind = grid.Face_Index(0, Vector3i(bx * 4 + (idx + 1), by * 4 + idy, bz * 4 + idz));
		if (idx != 3) face_x[face_ind] = shared_cell_data[(idz * 4 + idy) * 4 + (idx + 1)] - cell_data;

		// y-axis faces
		face_ind = grid.Face_Index(1, Vector3i(bx * 4 + idx, by * 4 + (idy + 1), bz * 4 + idz));
		if (idy != 3) face_y[face_ind] = shared_cell_data[(idz * 4 + (idy + 1)) * 4 + idx] - cell_data;

		// z-axis faces
		face_ind = grid.Face_Index(2, Vector3i(bx * 4 + idx, by * 4 + idy, bz * 4 + (idz + 1)));
		if (idz != 3) face_z[face_ind] = shared_cell_data[((idz + 1) * 4 + idy) * 4 + idx] - cell_data;

		// x-axis shared faces
		face_ind = grid.Face_Index(0, Vector3i(bx * 4 + idx, by * 4 + idy, bz * 4 + idz));
		if (idx == 0) atomicAdd(face_x + face_ind, cell_data);

		face_ind = grid.Face_Index(0, Vector3i(bx * 4 + (idx + 1), by * 4 + idy, bz * 4 + idz));
		if (idx == 3) atomicAdd(face_x + face_ind, -cell_data);

		// y-axis shared faces
		face_ind = grid.Face_Index(1, Vector3i(bx * 4 + idx, by * 4 + idy, bz * 4 + idz));
		if (idy == 0) atomicAdd(face_y + face_ind, cell_data);

		face_ind = grid.Face_Index(1, Vector3i(bx * 4 + idx, by * 4 + (idy + 1), bz * 4 + idz));
		if (idy == 3) atomicAdd(face_y + face_ind, -cell_data);

		// z-axis shared faces
		face_ind = grid.Face_Index(2, Vector3i(bx * 4 + idx, by * 4 + idy, bz * 4 + idz));
		if (idz == 0) atomicAdd(face_z + face_ind, cell_data);

		face_ind = grid.Face_Index(2, Vector3i(bx * 4 + idx, by * 4 + idy, bz * 4 + (idz + 1)));
		if (idz == 3) atomicAdd(face_z + face_ind, -cell_data);
	}

	//diff of C, store on F
	template<class T, int d>
	void D_CoCell_Mapping(const Field<T,d,DataHolder::DEVICE> &C, FaceField<T,d,DataHolder::DEVICE> &F)
	{
		F.Fill(0);
		if constexpr (d == 2) {
			int Nx = C.grid.counts[0], Ny = C.grid.counts[1];
			T* face_x = thrust::raw_pointer_cast(F.face_data[0].data());
			T* face_y = thrust::raw_pointer_cast(F.face_data[1].data());
			const T* cell = thrust::raw_pointer_cast(C.data.data());

			D_CoCell_Mapping_Kernel2 << <dim3((Nx >> 3), (Ny >> 3)), dim3(8, 8) >> > (C.grid, face_x, face_y, cell);
		}
		else if constexpr (d == 3) {
			int Nx = C.grid.counts[0], Ny = C.grid.counts[1], Nz = C.grid.counts[2];
			T* face_x = thrust::raw_pointer_cast(F.face_data[0].data());
			T* face_y = thrust::raw_pointer_cast(F.face_data[1].data());
			T* face_z = thrust::raw_pointer_cast(F.face_data[2].data());
			const T* cell = thrust::raw_pointer_cast(C.data.data());

			D_CoCell_Mapping_Kernel3 << <dim3((Nx >> 2), (Ny >> 2), (Nz >> 2)), dim3(4, 4, 4) >> > (C.grid, face_x, face_y, face_z, cell);
		}
	}

	// for blockDim = (8, 8)
	// iterate through cell
	// velocity in face_x, face_y
	// calculate divergence to cell
	template<class T>
	__global__ void D_Face_Mapping_Kernel2(Grid<2,CENTER> grid, T* cell, const T* face_x, const T* face_y)
	{
		const int nbx = gridDim.x;
		const int nby = gridDim.y;

		const int bx = blockIdx.x;
		const int by = blockIdx.y;

		const int idx = threadIdx.x;
		const int idy = threadIdx.y;

		__shared__ T shared_face_data[72];

		T div;
		{
			// 9x8 face_x data load
			shared_face_data[idy * 8 + idx] = face_x[grid.Face_Index(0, Vector2i(bx * 8 + idy, by * 8 + idx))];
			if (idy == 0) shared_face_data[64 + idx] = face_x[grid.Face_Index(0, Vector2i((bx + 1) * 8 + idy, by * 8 + idx))];
			__syncthreads();

			// left x-axis faces
			div = -shared_face_data[idx * 8 + idy];

			// right x-axis faces
			div += shared_face_data[(idx + 1) * 8 + idy];
		}
		__syncthreads();

		{
			// 8x9 face_y data load
			shared_face_data[idy * 8 + idx] = face_y[grid.Face_Index(1, Vector2i(bx * 8 + idx, by * 8 + idy))];
			if (idy == 0) shared_face_data[64 + idx] = face_y[grid.Face_Index(1, Vector2i(bx * 8 + idx, (by + 1) * 8 + idy))];
			__syncthreads();

			// down y-axis faces
			div -= shared_face_data[idy * 8 + idx];

			// up y-axis faces
			div += shared_face_data[(idy + 1) * 8 + idx];
		}

		cell[grid.Index(Vector2i(bx * 8 + idx, by * 8 + idy))] = div;
	}

	// for blockDim = (4, 4, 4)
	// iterate through cell
	template<class T>
	__global__ void D_Face_Mapping_Kernel3(Grid<3, CENTER> grid, T* cell, const T* face_x, const T* face_y, const T* face_z)
	{
		const int nbx = gridDim.x;
		const int nby = gridDim.y;
		const int nbz = gridDim.z;

		const int bx = blockIdx.x;
		const int by = blockIdx.y;
		const int bz = blockIdx.z;

		const int idx = threadIdx.x;
		const int idy = threadIdx.y;
		const int idz = threadIdx.z;

		__shared__ T shared_face_data[80];

		T div = 0;
		{
			// 5x4x4 face_x data load
			shared_face_data[(idz * 4 + idy) * 4 + idx] = face_x[grid.Face_Index(0, Vector3i(bx * 4 + idz, by * 4 + idx, bz * 4 + idy))];
			if (idz == 0) shared_face_data[(4 * 4 + idy) * 4 + idx] = face_x[grid.Face_Index(0, Vector3i(bx * 4 + 4, by * 4 + idx, bz * 4 + idy))];
			__syncthreads();

			// left x-axis faces
			div = -shared_face_data[(idx * 4 + idz) * 4 + idy];

			// right x-axis faces
			div += shared_face_data[((idx + 1) * 4 + idz) * 4 + idy];
		}
		__syncthreads();

		{
			// 4x5x4 face_y data load
			shared_face_data[(idz * 4 + idy) * 4 + idx] = face_y[grid.Face_Index(1, Vector3i(bx * 4 + idx, by * 4 + idz, bz * 4 + idy))];
			if (idz == 0) shared_face_data[(4 * 4 + idy) * 4 + idx] = face_y[grid.Face_Index(1, Vector3i(bx * 4 + idx, by * 4 + 4, bz * 4 + idy))];
			__syncthreads();

			// down y-axis faces
			div += -shared_face_data[(idy * 4 + idz) * 4 + idx];

			// up y-axis faces
			div += shared_face_data[((idy + 1) * 4 + idz) * 4 + idx];
		}
		__syncthreads();

		{
			// 4x4x5 face_z data load
			shared_face_data[(idz * 4 + idy) * 4 + idx] = face_z[grid.Face_Index(2, Vector3i(bx * 4 + idx, by * 4 + idy, bz * 4 + idz))];
			if (idz == 0) shared_face_data[(4 * 4 + idy) * 4 + idx] = face_z[grid.Face_Index(2, Vector3i(bx * 4 + idx, by * 4 + idy, bz * 4 + 4))];
			__syncthreads();

			// front z-axis faces
			div += -shared_face_data[(idz * 4 + idy) * 4 + idx];

			// back y-axis faces
			div += shared_face_data[((idz + 1) * 4 + idy) * 4 + idx];
		}
		__syncthreads();

		cell[grid.Index(Vector3i(bx * 4 + idx, by * 4 + idy, bz * 4 + idz))] = div;
	}

	template<class T, int d>
	void D_Face_Mapping(const FaceField<T, d, DEVICE>& F, Field<T, d, DEVICE>& C) {
		if constexpr (d == 2) {
			int Nx = F.grid.counts[0], Ny = F.grid.counts[1];
			T* cell = thrust::raw_pointer_cast(C.data.data());
			const T* face_x = thrust::raw_pointer_cast(F.face_data[0].data());
			const T* face_y = thrust::raw_pointer_cast(F.face_data[1].data());
			D_Face_Mapping_Kernel2 << <dim3((Nx >> 3), (Ny >> 3)), dim3(8, 8) >> > (F.grid, cell, face_x, face_y);
		}
		else if constexpr (d == 3) {
			int Nx = F.grid.counts[0], Ny = F.grid.counts[1], Nz = F.grid.counts[2];
			T* cell = thrust::raw_pointer_cast(C.data.data());
			const T* face_x = thrust::raw_pointer_cast(F.face_data[0].data());
			const T* face_y = thrust::raw_pointer_cast(F.face_data[1].data());
			const T* face_z = thrust::raw_pointer_cast(F.face_data[2].data());
			D_Face_Mapping_Kernel3 << <dim3((Nx >> 2), (Ny >> 2), (Nz >> 2)), dim3(4, 4, 4) >> > (F.grid, cell, face_x, face_y, face_z);
		}
	}
}