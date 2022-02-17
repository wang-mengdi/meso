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
	__global__ void D_CoCell_Mapping_Kernel(Grid<2> grid, T* face_x, T* face_y, const T* cell)
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

	//diff of C, store on F
	template<class T, int d>
	void D_CoCell_Mapping(const Field<T,d,DataHolder::DEVICE> &C, FaceField<T,d,DataHolder::DEVICE> &F)
	{
		if constexpr (d == 2) {
			int Nx = C.grid.counts[0], Ny = C.grid.counts[1];
			T* face_x = thrust::raw_pointer_cast(F.face_data[0].data());
			T* face_y = thrust::raw_pointer_cast(F.face_data[1].data());
			const T* cell = thrust::raw_pointer_cast(C.data.data());

			D_CoCell_Mapping_Kernel << <dim3((Nx >> 3), (Ny >> 3)), dim3(8, 8) >> > (C.grid, face_x, face_y, cell);
		}
	}

	// for blockDim = (8, 8)
	// iterate through cell
	// velocity in face_x, face_y
	// calculate divergence to cell
	template<class T>
	__global__ void D_Face_Mapping_Kernel(Grid<2,GridType::CELL> grid, T* cell, const T* face_x, const T* face_y)
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
			div += shared_face_data[idy * 8 + idx];

			// up y-axis faces
			div += -shared_face_data[(idy + 1) * 8 + idx];
		}

		cell[grid.Index(Vector2i(bx * 8 + idx, by * 8 + idy))] = div;
	}

	template<class T, int d>
	void D_Face_Mapping(const FaceField<T, d, DEVICE>& F, Field<T, d, DEVICE>& C) {
		if constexpr (d == 2) {
			int Nx = F.grid.counts[0], Ny = F.grid.counts[1];
			T* cell = thrust::raw_pointer_cast(C.data.data());
			const T* face_x = thrust::raw_pointer_cast(F.face_data[0].data());
			const T* face_y = thrust::raw_pointer_cast(F.face_data[1].data());
			D_Face_Mapping_Kernel << <dim3((Nx >> 3), (Ny >> 3)), dim3(8, 8) >> > (F.grid, cell, face_x, face_y);
		}
	}
}