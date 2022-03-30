//////////////////////////////////////////////////////////////////////////
// Basic grid data representation
// Copyright (c) (2022-), Bo Zhu, Zangyueyang Xian, Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "Common.h"
#include "AuxFunc.h"
#include "cuda_runtime.h"

namespace Meso {

	//Under CELL mode, positions are cell centers
	//Under NODE mode, positions are nodes
	enum GridType { NODE = 0, CELL };

	template<int d, GridType grid_type = GridType::CELL>
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

		__host__ __device__ VectorDi Counts(void) { return counts; }
		__host__ __device__ VectorDi Face_Counts(const int axis)const { VectorDi fcounts = counts; fcounts[axis] += block_size; return fcounts; }
		__host__ __device__ int DoF(void) const { return counts.prod(); }
		__host__ __device__ int Face_DoF(int axis)const { return Face_Counts(axis).prod(); }
		__host__ __device__ bool Valid(const VectorDi coord)const { bool res = true; for (int i = 0; i < d; i++) { res &= (0 <= coord[i] && coord[i] < counts[i]); }return res; }

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

		__host__ __device__ int Face_Index(const int axis, const VectorDi face) const
		{
			if constexpr (d == 2) {
				int x = face[0], y = face[1];
				int b_ind = (y >> 3) * ((counts[0] >> 3) + (axis == 0)) + (x >> 3);
				int idx = x & 7, idy = y & 7;
				return b_ind * 64 + ((axis == 0) * (idx * 8 + idy) + (axis == 1) * (idy * 8 + idx));
			}
			else if constexpr (d == 3) {
				int x = face[0], y = face[1], z = face[2];
				int nbx = (counts[0] >> 2) + (axis == 0), nby = (counts[1] >> 2) + (axis == 1), nbz = (counts[2] >> 2) + (axis == 2);
				int bx = x >> 2, by = y >> 2, bz = z >> 2;
				int idx = x & 0b11, idy = y & 0b11, idz = z & 0b11;
				return ((bz * nby + by) * nbx + bx) * 64 + (axis == 0) * ((idx * 4 + idz) * 4 + idy) + (axis == 1) * ((idy * 4 + idz) * 4 + idx) + (axis == 2) * ((idz * 4 + idy) * 4 + idx);
			}
		}

		__host__ __device__ VectorDi Face_Coord(const int axis, const int face_index)const {
			VectorDi face;
			if constexpr (d == 2) {
				int b_ind = face_index >> 6;
				int i1 = (face_index & 0b111000) >> 3, i0 = face_index & 0b111;
				int nbx = (counts[0] >> 3) + (axis == 0);
				int idx, idy;
				if (axis == 0) idx = i1, idy = i0;
				else idx = i0, idy = i1;
				face[0] = ((b_ind % nbx) << 3) + idx;
				face[1] = ((b_ind / nbx) << 3) + idy;
			}
			else if constexpr (d == 3) {
				int nbx = (counts[0] >> 2) + (axis == 0), nby = (counts[1] >> 2) + (axis == 1), nbz = (counts[2] >> 2) + (axis == 2);
				int idx = face_index & 0b11, idy = (face_index & 0b1100) >> 2, idz = (face_index & 0b110000) >> 4;
				int b = face_index >> 6;

				face[0] = ((b % nbx) << 2) + (axis == 0) * idz + (axis == 1) * idx + (axis == 2) * idx;
				face[1] = (((b / nbx) % nby) << 2) + (axis == 0) * idx + (axis == 1) * idz + (axis == 2) * idy;
				face[2] = ((b / nbx / nby) << 2) + (axis == 0) * idy + (axis == 1) * idy + (axis == 2) * idz;
			}
			return face;
		}

		////parallel iterators
		template<class Fcell>//Fcell is a (void) function takes a cell index
		void Exec_Cells(Fcell f) const {
			const int dof = DoF();
#pragma omp parallel for
			for (int c = 0; c < dof; c++) {
				const VectorDi cell = Coord(c);
				f(cell);
			}
		}
	};

}