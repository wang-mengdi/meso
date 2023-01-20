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


	};

}
