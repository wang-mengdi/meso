//////////////////////////////////////////////////////////////////////////
// Basic grid data representation
// Copyright (c) (2022-), Bo Zhu, Zangyueyang Xian, Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "Common.h"
#include "AuxFunc.h"
#include "cuda_runtime.h"

//Under CELL mode, positions are cell centers
//Under NODE mode, positions are nodes
enum GridType { NODE = 0, CELL};

template<int d, GridType grid_type=GridType::CELL>
class Grid {
	Typedef_VectorD(d);
public:
	static constexpr int block_size = (d == 2 ? 8 : 4);
	VectorDi counts;
	real dx;
	VectorD pos_min;

	Grid(const VectorDi _counts = VectorDi::Zero(), const real _dx = 0, const VectorD domain_min = VectorD::Zero()) :
		dx(_dx)
	{
		counts = VectorFunc::Round_Up_To_Align<d>(_counts, block_size);
		if (counts != _counts) Warn("Grid size not divisible by {} in dimension {}, automtically round up to {}", block_size, d, counts);
		if (grid_type == GridType::CELL) pos_min = domain_min;
		else pos_min = domain_min + VectorFunc::V<d>(0.5, 0.5, 0.5) * dx;
	}

	__host__ __device__ int DoF(void) const { return counts.prod(); }

	__host__ __device__ int Index(const VectorDi coord) const {
		if constexpr (d == 2) {
			return ((coord[1] >> 3) * (counts[0] >> 3) + (coord[0] >> 3)) * 64 + ((coord[1] & 7) * 8 + (coord[0] & 7));
		}
		else if constexpr (d == 3) {
			int nbx = counts[0] >> 2, nby = counts[1] >> 2, nbz = counts[2] >> 2;
			int bx = coord[0] >> 2, by = coord[1] >> 2, bz = coord[2] >> 2;
			int idx = coord[0] & 0b11, idy = coord[1] & 0b11, idz = coord[2] & 0b11;
			return ((bz * nby + by) * nbx + bx) * 64 + ((idz * 4 + idy) * 4 + idx);
		}
		else Assert(false, "Grid::Index: d==2 or d==3");
	}

	__host__ __device__ VectorDi Coord(const int index)const {
		if constexpr (d == 2) {
			int idx = index & 0b111;
			int idy = (index & 0b111000) >> 3;
			int b = index >> 6;
			int x = ((b % (counts[0] >> 3)) << 3) + idx;
			int y = ((b / (counts[0] >> 3)) << 3) + idy;
			return Vector2i(x, y);
		}
		else if constexpr (d == 3) {
			int nbx = counts[0] >> 2, nby = counts[1] >> 2, nbz = counts[2] >> 2;
			int idx = index & 0b11, idy = (index & 0b1100) >> 2, idz = (index & 0b110000) >> 4;
			int b = index >> 6;
			int x = ((b % nbx) << 2) + idx;
			int y = (((b / nbx) % nby) << 2) + idy;
			int z = ((b / nbx / nby) << 2) + idz;
			return Vector3i(x, y, z);
		}
	}

	////parallel iterators
	template<class Fcell>//Fcell is a (void) function takes a cell index
	void Exec_Each(Fcell f) const {
		const int dof = DoF();
#pragma omp parallel for
		for (int c = 0; c < dof; c++) {
			const VectorDi cell = Coord(c);
			f(cell);
		}
	}
};