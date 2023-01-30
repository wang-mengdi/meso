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
	__global__ void Poisson_Boundary_Apply_Kernel(const Grid<d> _grid, const int* _boundary_tiles, unsigned char* _cell_type, T** _vol, const T* _p, T* _Ap)
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
	__global__ void Poisson_Apply_Kernel2(const Grid<2> _grid, const unsigned char* _cell_type, T** _vol,
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
			if (shared_cell_type[idx][idy + 1] == 1 || shared_cell_type[idx][idy + 1] == 2)
				result += cell_p * shared_vol_x[idx][idy];
			else
				result += (cell_p - shared_p[idx][idy + 1]) * shared_vol_x[idx][idy];
			if (shared_cell_type[idx + 2][idy + 1] == 1 || shared_cell_type[idx + 2][idy + 1] == 2)
				result += cell_p * shared_vol_x[idx + 1][idy];
			else
				result += (cell_p - shared_p[idx + 2][idy + 1]) * shared_vol_x[idx + 1][idy];
			if (shared_cell_type[idx + 1][idy] == 1 || shared_cell_type[idx + 1][idy] == 2)
				result += cell_p * shared_vol_y[idx][idy];
			else
				result += (cell_p - shared_p[idx + 1][idy]) * shared_vol_y[idx][idy];
			if (shared_cell_type[idx + 1][idy + 2] == 1 || shared_cell_type[idx + 1][idy + 2] == 2)
				result += cell_p * shared_vol_y[idx][idy + 1];
			else
				result += (cell_p - shared_p[idx + 1][idy + 2]) * shared_vol_y[idx][idy + 1];
		}
		_Ap[global_id] = result;
	}

	template<class T>
	__global__ void Poisson_Apply_Kernel3(const Grid<3> _grid, const unsigned char* _cell_type, T** _vol,
		const T* _p, T* _Ap)
	{
		Typedef_VectorD(3);
		// calculate index
		VectorDi coord = GPUFunc::Thread_Coord<3>(blockIdx, threadIdx);
		const int idx = threadIdx.x;
		const int idy = threadIdx.y;
		const int idz = threadIdx.z;
		const int global_id = _grid.Index(coord);
		const int grid_counts_x = _grid.Counts()[0];
		const int grid_counts_y = _grid.Counts()[1];
		const int grid_counts_z = _grid.Counts()[2];

		// define shared memory
		__shared__ unsigned char shared_cell_type[6][6][6];
		__shared__ T shared_p[6][6][6];
		__shared__ T shared_vol_x[5][4][4];
		__shared__ T shared_vol_y[4][5][4];
		__shared__ T shared_vol_z[4][4][5];
		
		// load data to shared memory
		unsigned char type = _cell_type[global_id];
		shared_cell_type[idx + 1][idy + 1][idz + 1] = type;
		if (idx == 0)
		{
			if (coord[0] == 0)
				shared_cell_type[0][idy + 1][idz + 1] = 1;
			else
				shared_cell_type[0][idy + 1][idz + 1] = _cell_type[_grid.Index(coord - VectorDi::Unit(0))];
		}
		if (idx == 3)
		{
			if (coord[0] == grid_counts_x - 1)
				shared_cell_type[5][idy + 1][idz + 1] = 1;
			else
				shared_cell_type[5][idy + 1][idz + 1] = _cell_type[_grid.Index(coord + VectorDi::Unit(0))];
		}
		if (idy == 0)
		{
			if (coord[1] == 0)
				shared_cell_type[idx + 1][0][idz + 1] = 1;
			else
				shared_cell_type[idx + 1][0][idz + 1] = _cell_type[_grid.Index(coord - VectorDi::Unit(1))];
		}
		if (idy == 3)
		{
			if (coord[1] == grid_counts_y - 1)
				shared_cell_type[idx + 1][5][idz + 1] = 1;
			else
				shared_cell_type[idx + 1][5][idz + 1] = _cell_type[_grid.Index(coord + VectorDi::Unit(1))];
		}
		if (idz == 0)
		{
			if (coord[2] == 0)
				shared_cell_type[idx + 1][idy + 1][0] = 1;
			else
				shared_cell_type[idx + 1][idy + 1][0] = _cell_type[_grid.Index(coord - VectorDi::Unit(2))];
		}
		if (idz == 3)
		{
			if (coord[2] == grid_counts_z - 1)
				shared_cell_type[idx + 1][idy + 1][5] = 1;
			else
				shared_cell_type[idx + 1][idy + 1][5] = _cell_type[_grid.Index(coord + VectorDi::Unit(2))];
		}
		
		T cell_p = _p[global_id];
		shared_p[idx + 1][idy + 1][idz + 1] = cell_p;
		if (idx == 0)
		{
			if (coord[0] == 0)
				shared_p[0][idy + 1][idz + 1] = 0;
			else
				shared_p[0][idy + 1][idz + 1] = _p[_grid.Index(coord - VectorDi::Unit(0))];
		}
		if (idx == 3)
		{
			if (coord[0] == grid_counts_x - 1)
				shared_p[5][idy + 1][idz + 1] = 0;
			else
				shared_p[5][idy + 1][idz + 1] = _p[_grid.Index(coord + VectorDi::Unit(0))];
		}
		if (idy == 0)
		{
			if (coord[1] == 0)
				shared_p[idx + 1][0][idz + 1] = 0;
			else
				shared_p[idx + 1][0][idz + 1] = _p[_grid.Index(coord - VectorDi::Unit(1))];
		}
		if (idy == 3)
		{
			if (coord[1] == grid_counts_y - 1)
				shared_p[idx + 1][5][idz + 1] = 0;
			else
				shared_p[idx + 1][5][idz + 1] = _p[_grid.Index(coord + VectorDi::Unit(1))];
		}
		if (idz == 0)
		{
			if (coord[2] == 0)
				shared_p[idx + 1][idy + 1][0] = 0;
			else
				shared_p[idx + 1][idy + 1][0] = _p[_grid.Index(coord - VectorDi::Unit(2))];
		}
		if (idz == 3)
		{
			if (coord[2] == grid_counts_z - 1)
				shared_p[idx + 1][idy + 1][5] = 0;
			else
				shared_p[idx + 1][idy + 1][5] = _p[_grid.Index(coord + VectorDi::Unit(2))];
		}
		
		shared_vol_x[idx][idy][idz] = _vol[0][_grid.Face_Index(0, coord)];
		if (idx == 3)
			shared_vol_x[4][idy][idz] = _vol[0][_grid.Face_Index(0, coord + VectorDi::Unit(0))];

		shared_vol_y[idx][idy][idz] = _vol[1][_grid.Face_Index(1, coord)];
		if (idy == 3)
			shared_vol_y[idx][4][idz] = _vol[1][_grid.Face_Index(1, coord + VectorDi::Unit(1))];
		
		shared_vol_z[idx][idy][idz] = _vol[2][_grid.Face_Index(2, coord)];
		if (idz == 3)
			shared_vol_z[idx][idy][4] = _vol[2][_grid.Face_Index(2, coord + VectorDi::Unit(2))];

		__syncthreads();
		
		T result = 0;
		if (type == 1 || type == 2)
			result = cell_p;
		else
		{
			if (shared_cell_type[idx][idy + 1][idz + 1] == 1 || shared_cell_type[idx][idy + 1][idz + 1] == 2)
				result += cell_p * shared_vol_x[idx][idy][idz];
			else
				result += (cell_p - shared_p[idx][idy + 1][idz + 1]) * shared_vol_x[idx][idy][idz];
			if (shared_cell_type[idx + 2][idy + 1][idz + 1] == 1 || shared_cell_type[idx + 2][idy + 1][idz + 1] == 2)
				result += cell_p * shared_vol_x[idx + 1][idy][idz];
			else
				result += (cell_p - shared_p[idx + 2][idy + 1][idz + 1]) * shared_vol_x[idx + 1][idy][idz];
			if (shared_cell_type[idx + 1][idy][idz + 1] == 1 || shared_cell_type[idx + 1][idy][idz + 1] == 2)
				result += cell_p * shared_vol_y[idx][idy][idz];
			else
				result += (cell_p - shared_p[idx + 1][idy][idz + 1]) * shared_vol_y[idx][idy][idz];
			if (shared_cell_type[idx + 1][idy + 2][idz + 1] == 1 || shared_cell_type[idx + 1][idy + 2][idz + 1] == 2)
				result += cell_p * shared_vol_y[idx][idy + 1][idz];
			else
				result += (cell_p - shared_p[idx + 1][idy + 2][idz + 1]) * shared_vol_y[idx][idy + 1][idz];
			if (shared_cell_type[idx + 1][idy + 1][idz] == 1 || shared_cell_type[idx + 1][idy + 1][idz] == 2)
				result += cell_p * shared_vol_z[idx][idy][idz];
			else
				result += (cell_p - shared_p[idx + 1][idy + 1][idz]) * shared_vol_z[idx][idy][idz];
			if (shared_cell_type[idx + 1][idy + 1][idz + 2] == 1 || shared_cell_type[idx + 1][idy + 1][idz + 2] == 2)
				result += cell_p * shared_vol_z[idx][idy][idz + 1];
			else
				result += (cell_p - shared_p[idx + 1][idy + 1][idz + 2]) * shared_vol_z[idx][idy][idz + 1];
		}
		_Ap[global_id] = result;
	}

	//Negative Poisson mapping -lap(p), except some masked points
	//Masked cells will be viewed as 0 in poisson mapping
	//Which means adjacent faces of masked cells will have volume 0
	template<class T, int d>
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

		MaskedPoissonMapping() {}
		MaskedPoissonMapping(const Grid<d> grid) { Allocate_Memory(grid); }

		void Allocate_Memory(const Grid<d> grid) {
			Assert(grid.Is_Unpadded(), "MaskedPoissonMapping: invalid grid {}, padding not allowed", grid);
			dof = grid.Memory_Size();
			cell_type.Init(grid);
			vol.Init(grid);
			is_boundary_tile.resize(grid.Memory_Size() / 64);
		}

		void Init(const Grid<d> grid) {
			Allocate_Memory(grid);
		}

		template<DataHolder side>
		void Init(const Field<unsigned char, d, side>& _cell_type, const FaceField<T, d, side>& _vol) {
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
			if constexpr (d == 2)
				vol.grid.Exec_Kernel(Poisson_Apply_Kernel2<T>, vol.grid, cell_type.Data_Ptr(), ArrayFunc::Data(vol.face_data_ptr),
					ArrayFunc::Data(p), ArrayFunc::Data(Ap));
			else
				vol.grid.Exec_Kernel(Poisson_Apply_Kernel3<T>, vol.grid, cell_type.Data_Ptr(), ArrayFunc::Data(vol.face_data_ptr),
					ArrayFunc::Data(p), ArrayFunc::Data(Ap));
		}

		void Boundary_Apply(ArrayDv<T>& Ap, const ArrayDv<T>& p)
		{
			Poisson_Boundary_Apply_Kernel<T, d> << < boundary_tiles.size(), 64 >> > (vol.grid, ArrayFunc::Data(boundary_tiles),
				cell_type.Data_Ptr(), ArrayFunc::Data(vol.face_data_ptr), ArrayFunc::Data(p), ArrayFunc::Data(Ap));
		}
	};

}
