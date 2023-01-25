//////////////////////////////////////////////////////////////////////////
// Poisson Linear Mapping with celltype and vol
// Copyright (c) (2022-), Zangyueyang Xian, Mengdi Wang, Yuchen Sun
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "LinearMapping.h"
#include "Field.h"
#include "FaceField.h"
#include "ExteriorDerivative.h"
#include "AuxFunc.h"
#include "BoundaryCondition.h"
using namespace thrust::placeholders;
namespace Meso {
	template<int d>
	__global__ void Refresh_Boundary(const unsigned char _cur_dis, const Grid<d> _grid, unsigned char* _cell_type_ptr)
	{
		Typedef_VectorD(d);
		VectorDi coord = GPUFunc::Thread_Coord<d>(blockIdx, threadIdx);
		int idx = _grid.Index(coord);
		unsigned char type = _cell_type_ptr[idx] & 0b11;
		if (type != 0)
			return;
		if (_cur_dis == 0)
			for (int axis = 0; axis < d; axis++)
				for (int side = -1; side <= 1; side += 2)
				{
					VectorDi  nb_coord = coord + VectorDi::Unit(axis) * side;
					if (_grid.Valid(nb_coord))
					{
						int nb_idx = _grid.Index(nb_coord);
						unsigned char nb_type = _cell_type_ptr[nb_idx] & 0b11;
						if (nb_type == 1 || nb_type == 2)
						{
							_cell_type_ptr[idx] = 3;
							return;
						}
					}
					else
					{
						_cell_type_ptr[idx] = 3;
						return;
					}
				}
		else
			for (int axis = 0; axis < d; axis++)
				for (int side = -1; side <= 1; side += 2)
				{
					VectorDi  nb_coord = coord + VectorDi::Unit(axis) * side;
					if (_grid.Valid(nb_coord))
					{
						int nb_idx = _grid.Index(nb_coord);
						unsigned nb_flag = _cell_type_ptr[nb_idx];
						unsigned char nb_type = nb_flag & 0b11;
						if (nb_type == 3)
						{
							unsigned char nb_dis = nb_flag >> 2;
							if (nb_dis == _cur_dis - 1)
							{
								_cell_type_ptr[idx] = (_cur_dis << 2) + 0b11;
								return;
							}
						}
					}
				}
	}

	template<int d>
	__global__ void Clear_Dis(const Grid<d> _grid, unsigned char* _cell_type_ptr)
	{
		Typedef_VectorD(d);
		VectorDi coord = GPUFunc::Thread_Coord<d>(blockIdx, threadIdx);
		int idx = _grid.Index(coord);
		_cell_type_ptr[idx] = _cell_type_ptr[idx] & 0b11;
	}

	template<int d>
	__global__ void Mark_Boundary_Tile(const Grid<d> _grid, const unsigned char* _cell_type, int* _is_boundary_tile)
	{
		Typedef_VectorD(d);
		VectorDi coord = GPUFunc::Thread_Coord<d>(blockIdx, threadIdx);
		int idx = _grid.Index(coord);
		if (_cell_type[idx] == 3)
		{
			int block_idx = idx >> 6;
			_is_boundary_tile[block_idx] = 1;
		}
	}

	template<class T, int d>
	__global__ void Poisson_Boundary_Apply_Kernel(const Grid<d> _grid, const int* _boundary_tiles, unsigned char* _cell_type, const T** _vol, const T* _p, T* _Ap)
	{
		Typedef_VectorD(d);
		int block_id = _boundary_tiles[blockIdx.x];
		int cell_id = block_id * blockDim.x + threadIdx.x;
		if (_cell_type[cell_id] != 3)
			return;
		VectorDi cell = _grid.Coord(cell_id);
		T cell_p = _p[cell_id];
		T cell_Ap = 0;
		for (int axis = 0; axis < d; axis++)
			for (int side = -1; side <= 1; side += 2)
			{
				VectorDi nb_cell = cell + VectorDi::Unit(axis) * side;
				int nb_cell_id = _grid.Index(nb_cell);
				T nb_cell_p;
				if (!_grid.Valid(nb_cell) || _cell_type[nb_cell_id] == 1 || _cell_type[nb_cell_id] == 2)
					nb_cell_p = 0;
				else
					nb_cell_p = _p[nb_cell_id];
				T face_vol;
				if (side == -1)
					face_vol = _vol[axis][_grid.Face_Index(axis, cell)];
				else
					face_vol = _vol[axis][_grid.Face_Index(axis, nb_cell)];
				cell_Ap += face_vol * (cell_p - nb_cell_p);
			}
		_Ap[cell_id] = cell_Ap;
	}

	template<class T>
	__global__ void Poisson_Apply_Kernel2(const Grid<2> _grid, const unsigned char* _cell_type, const T** _vol,
		const T* _p, T* _Ap)
	{
		Typedef_VectorD(2);
		// calculate index
		VectorDi coord = GPUFunc::Thread_Coord<2>(blockIdx, threadIdx);
		const int idx = threadIdx.x;
		const int idy = threadIdx.y;
		const int global_id = _grid.Index(coord);
		const int grid_counts_x = _grid.Counts()[0];
		const int grid_counts_y = _grid.Counts()[1];

		// define shared memory
		__shared__ unsigned char shared_cell_type[10][10];
		__shared__ T shared_p[10][10];
		__shared__ T shared_vol_x[9][8];
		__shared__ T shared_vol_y[8][9];

		// load data to shared memory
		unsigned char type = _cell_type[global_id];
		shared_cell_type[idx + 1][idy + 1] = type;
		if (idx == 0)
		{
			if (coord[0] == 0)
				shared_cell_type[0][idy + 1] = 1;
			else
				shared_cell_type[0][idy + 1] = _cell_type[_grid.Index(coord-VectorDi::Unit(0))];
		}
		if (idx == 7)
		{
			if (coord[0] == grid_counts_x - 1)
				shared_cell_type[9][idy + 1] = 1;
			else
				shared_cell_type[9][idy + 1] = _cell_type[_grid.Index(coord + VectorDi::Unit(0))];
		}
		if (idy == 0)
		{
			if (coord[1] == 0)
				shared_cell_type[idx + 1][0] = 1;
			else
				shared_cell_type[idx + 1][0] = _cell_type[_grid.Index(coord - VectorDi::Unit(1))];
		}
		if (idy == 7)
		{
			if (coord[1] == grid_counts_y - 1)
				shared_cell_type[idx + 1][9] = 1;
			else
				shared_cell_type[idx + 1][9] = _cell_type[_grid.Index(coord + VectorDi::Unit(1))];
		}
		__syncthreads();

		T cell_p = _p[global_id];
		shared_p[idx + 1][idy + 1] = cell_p;
		if (idx == 0)
		{
			if (coord[0] == 0)
				shared_p[0][idy + 1] = 0;
			else
				shared_p[0][idy + 1] = _p[_grid.Index(coord - VectorDi::Unit(0))];
		}
		if (idx == 7)
		{
			if (coord[0] == grid_counts_x - 1)
				shared_p[9][idy + 1] = 0;
			else
				shared_p[9][idy + 1] = _p[_grid.Index(coord + VectorDi::Unit(0))];
		}
		if (idy == 0)
		{
			if (coord[1] == 0)
				shared_p[idx + 1][0] = 0;
			else
				shared_p[idx + 1][0] = _p[_grid.Index(coord - VectorDi::Unit(1))];
		}
		if (idy == 7)
		{
			if (coord[1] == grid_counts_y - 1)
				shared_p[idx + 1][9] = 0;
			else
				shared_p[idx + 1][9] = _p[_grid.Index(coord + VectorDi::Unit(1))];
		}

		if (type == 1 || type == 2)
			shared_p[idx + 1][idy + 1] = 0;
		if (idx == 0)
		{
			unsigned nb_type = shared_cell_type[idx][idy + 1];
			if (nb_type == 1 || nb_type == 2)
				shared_p[idx][idy + 1] = 0;
		}
		if (idx == 7)
		{
			unsigned nb_type = shared_cell_type[idx + 2][idy + 1];
			if (nb_type == 1 || nb_type == 2)
				shared_p[idx + 2][idy + 1] = 0;
		}
		if (idy == 0)
		{
			unsigned nb_type = shared_cell_type[idx + 1][idy];
			if (nb_type == 1 || nb_type == 2)
				shared_p[idx + 1][idy] = 0;
		}
		if (idy == 7)
		{
			unsigned nb_type = shared_cell_type[idx + 1][idy + 2];
			if (nb_type == 1 || nb_type == 2)
				shared_p[idx + 1][idy + 2] = 0;
		}

		shared_vol_x[idx][idy] = _vol[0][_grid.Face_Index(0, coord)];
		if (idx == 7)
			shared_vol_x[8][idy] = _vol[0][_grid.Face_Index(0, coord + VectorDi::Unit(0))];

		shared_vol_y[idx][idy] = _vol[1][_grid.Face_Index(1, coord)];
		if (idy == 7)
			shared_vol_y[idx][8] = _vol[1][_grid.Face_Index(1, coord + VectorDi::Unit(1))];

		__syncthreads();

		T result = 0;
		if (type == 1 || type == 2)
			result = cell_p;
		else
		{
			result += (cell_p - shared_p[idx][idy + 1]) * shared_vol_x[idx][idy];
			result += (cell_p - shared_p[idx + 2][idy + 1]) * shared_vol_x[idx + 1][idy];
			result += (cell_p - shared_p[idx + 1][idy]) * shared_vol_y[idx][idy];
			result += (cell_p - shared_p[idx + 1][idy + 2]) * shared_vol_y[idx][idy + 1];
		}
		_Ap[global_id] = result;
	}

	//Negative Poisson mapping -lap(p), except some masked points
	//Masked cells will be viewed as 0 in poisson mapping
	//Which means adjacent faces of masked cells will have volume 0
	template<class T, int d, class ExteriorDerivative=ExteriorDerivativePadding0>
	class MaskedPoissonMapping : public LinearMapping<T> {
		using Base=LinearMapping<T>;
		//Ap=-lap(p)
		Typedef_VectorD(d);
	public:
		int dof;
		FieldDv<unsigned char, d> cell_type;
		FaceFieldDv<T, d> vol;
		ArrayDv<int> is_boundary_tile;
		ArrayDv<int> boundary_tiles;

		FaceFieldDv<T, d> temp_face;
		
		MaskedPoissonMapping() {}
		MaskedPoissonMapping(const Grid<d> grid) { Allocate_Memory(grid); }

		void Allocate_Memory(const Grid<d> grid) {
			Assert(grid.Is_Unpadded(), "MaskedPoissonMapping: invalid grid {}, padding not allowed", grid);
			dof = grid.Memory_Size();
			cell_type.Init(grid);
			vol.Init(grid);
			is_boundary_tile.resize(grid.Memory_Size() / 64);
			temp_face.Init(grid);
		}

		void Init(const Grid<d> grid) {
			Allocate_Memory(grid);
		}

		void Init(const Field<unsigned char, d>& _cell_type, const FaceField<T, d>& _vol) {
			Assert(_cell_type.grid.Indexer() == _vol.grid.Indexer(), "MaskedPoissonMapping::Init error: _cell_type grid {} unequal to _vol grid {}", _cell_type.grid, _vol.grid);
			Allocate_Memory(_cell_type.grid);
			cell_type.Deep_Copy(_cell_type);
			vol.Deep_Copy(_vol);

		}

		void Search_Boundary(const int _boundary_width = 3)
		{
			Meso::Grid<d> grid(cell_type.grid);
			for (int i = 0; i < _boundary_width; i++)
				grid.Exec_Kernel(Refresh_Boundary<d>, i, grid, cell_type.Data_Ptr());
			grid.Exec_Kernel(Clear_Dis<d>, grid, cell_type.Data_Ptr());
			ArrayFunc::Fill(is_boundary_tile, 0);
			grid.Exec_Kernel(Mark_Boundary_Tile<d>, grid, cell_type.Data_Ptr(), ArrayFunc::Data(is_boundary_tile));
			int num_boundary_tile = ArrayFunc::Sum<int>(is_boundary_tile);
			boundary_tiles.resize(num_boundary_tile);
			thrust::copy_if(thrust::make_counting_iterator<int>(0), thrust::make_counting_iterator<int>(grid.Memory_Size() / 64), 
				is_boundary_tile.begin(), boundary_tiles.begin(), thrust::identity<int>());
			Assert(boundary_tiles.size(), "boundary_tiles empty");
		}

		const Grid<d>& Grid(void) const {
			return vol.grid;
		}

		virtual int XDoF() const { return dof; }//number of cols

		virtual int YDoF() const { return dof; }//number of rows

		//input p, get Ap
		virtual void Apply(ArrayDv<T>& Ap, const ArrayDv<T>& p) {
			Base::Memory_Check(Ap, p, "PoissonMapping::Apply error: not enough space");

			auto fix = [=] __device__(T & tv,  unsigned char type) { if (type == 1 || type == 2) tv = (T)0; };
			auto multi = [=] __device__(T & tv, T tvol) { tv *= tvol; };
			auto neg = [=]__device__(T & tv) { tv = -tv; };
			auto cond_set = [=]__device__(T & tv1, T tv2, unsigned char type) { if (type == 1 || type == 2) tv1 = tv2; };

			Meso::Grid<d> grid = Grid();
			T* Ap_ptr = thrust::raw_pointer_cast(Ap.data());
			const T* p_ptr = thrust::raw_pointer_cast(p.data());

			//temp_cell=p, set to 0 for dirichlet and neumann cell
			cudaMemcpy(Ap_ptr, p_ptr, sizeof(T)* dof, cudaMemcpyDeviceToDevice);
			GPUFunc::Cwise_Mapping_Wrapper(Ap_ptr, cell_type.Data_Ptr(), fix, dof);

			//temp_face = grad(temp_cell) *. vol
			//d(p) ----- 1-form
			for (int axis = 0; axis < d; axis++)
				cudaMemset(temp_face.Data_Ptr(axis), 0, sizeof(T) * grid.Face_Grid(axis).Memory_Size());
			if constexpr (d == 2) grid.Exec_Kernel(&D_CoCell_Kernel2_Padding0<T>, grid, temp_face.Data_Ptr(0), temp_face.Data_Ptr(1), Ap_ptr);
			else if constexpr (d == 3) grid.Exec_Kernel(&D_CoCell_Kernel3_Padding0<T>, grid, temp_face.Data_Ptr(0), temp_face.Data_Ptr(1), temp_face.Data_Ptr(2), Ap_ptr);
			
			//d(p) *. vol ----- 1-form
			for (int axis = 0; axis < d; axis++)
				GPUFunc::Cwise_Mapping_Wrapper(temp_face.Data_Ptr(axis), vol.Data_Ptr(axis), multi, grid.Face_Grid(axis).Memory_Size());
			
			//div(temp_face)
			//d*d(p) *. vol ----- 3-form
			if constexpr (d == 2) grid.Exec_Kernel(D_Face_Kernel2_Padding0<T>, grid, Ap_ptr, temp_face.Data_Ptr(0), temp_face.Data_Ptr(1));
			else if constexpr (d == 3)grid.Exec_Kernel(D_Face_Kernel3_Padding0<T>, grid, Ap_ptr, temp_face.Data_Ptr(0), temp_face.Data_Ptr(1), temp_face.Data_Ptr(2));
			
			//neg div
			GPUFunc::Cwise_Mapping_Wrapper(Ap_ptr, neg, dof);
			
			//transfer dirichlet and neumann cell data to Ap
			GPUFunc::Cwise_Mapping_Wrapper(Ap_ptr, p_ptr, cell_type.Data_Ptr(), cond_set, dof);
		}

		void Apply_Kernel2(ArrayDv<T>& Ap, const ArrayDv<T>& p) {
			if constexpr (d == 2)
			{
				ArrayDv<const T*> vol_ptr(d);
				for (int axis = 0; axis < d; axis++)
					vol_ptr[axis] = vol.Data_Ptr(axis);
				vol.grid.Exec_Kernel(Poisson_Apply_Kernel2<T>, vol.grid, cell_type.Data_Ptr(), ArrayFunc::Data(vol_ptr),
					ArrayFunc::Data(p), ArrayFunc::Data(Ap));
			}
		}

		void Boundary_Apply(ArrayDv<T>& Ap, const ArrayDv<T>& p)
		{
			ArrayDv<const T*> vol_ptr(d);
			for (int axis = 0; axis < d; axis++)
				vol_ptr[axis] = vol.Data_Ptr(axis);
			Poisson_Boundary_Apply_Kernel<T, d> << < boundary_tiles.size(), 64 >> > (vol.grid, ArrayFunc::Data(boundary_tiles),
				cell_type.Data_Ptr(), ArrayFunc::Data(vol_ptr), ArrayFunc::Data(p), ArrayFunc::Data(Ap));
			checkCudaErrors(cudaGetLastError());
		}
	};

}
