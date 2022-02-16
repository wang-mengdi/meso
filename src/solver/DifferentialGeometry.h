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
	__global__ void Diff_Cell_Kernel(Grid<2> grid, T* face_x, T* face_y, T* cell)
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
		shared_cell_data[idy * 8 + idx] = cell[grid.cell_ind(bx * 8 + idx, by * 8 + idy)];
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

	//diff of C, store on F
	template<class T, int d>
	void Diff_Cell(const Field<T,d,DataHolder::DEVICE> &C, FaceField<T,d,DataHolder::DEVICE> &F)
	{
		if constexpr (d == 2) {
			int Nx = C.grid.counts[0], Ny = C.grid.counts[1];
			T* face_x = thrust::raw_pointer_cast(F.face_data[0].data());
			T* face_y = thrust::raw_pointer_cast(F.face_data[1].data());
			T* cell = thrust::raw_pointer_cast(C.data.data());

			Diff_Cell_Kernel << <dim3((Nx >> 3), (Ny >> 3)), dim3(8, 8) >> > (C.grid, face_x, face_y, cell);
		}
	}

	// for blockDim = (8, 8)
	// iterate through cell
	// velocity in face_x, face_y
	// calculate divergence to cell
	__global__ void D1MappingKernel(grid2D grid, Scalar* cell, Scalar* face_x, Scalar* face_y)
	{
		const int nbx = gridDim.x;
		const int nby = gridDim.y;

		const int bx = blockIdx.x;
		const int by = blockIdx.y;

		const int idx = threadIdx.x;
		const int idy = threadIdx.y;

		__shared__ Scalar shared_face_data[72];

		Scalar div;
		{
			// 9x8 face_x data load
			shared_face_data[idy * 8 + idx] = face_x[grid.face_ind(bx * 8 + idy, by * 8 + idx, 0)];
			if (idy == 0) shared_face_data[64 + idx] = face_x[grid.face_ind((bx + 1) * 8 + idy, by * 8 + idx, 0)];
			__syncthreads();

			// left x-axis faces
			div = -shared_face_data[idx * 8 + idy];

			// right x-axis faces
			div += shared_face_data[(idx + 1) * 8 + idy];
		}
		__syncthreads();

		{
			// 8x9 face_y data load
			shared_face_data[idy * 8 + idx] = face_y[grid.face_ind(bx * 8 + idx, by * 8 + idy, 1)];
			if (idy == 0) shared_face_data[64 + idx] = face_y[grid.face_ind(bx * 8 + idx, (by + 1) * 8 + idy, 1)];
			__syncthreads();

			// down y-axis faces
			div += shared_face_data[idy * 8 + idx];

			// up y-axis faces
			div += -shared_face_data[(idy + 1) * 8 + idx];
		}

		cell[grid.cell_ind(bx * 8 + idx, by * 8 + idy)] = div;
	}
}