//////////////////////////////////////////////////////////////////////////
// Basic grid data representation
// Copyright (c) (2022-), Bo Zhu, Zangyueyang Xian, Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "Common.h"
#include "AuxFunc.h"
//#include "cuda_runtime.h"

namespace Meso {

	//How to explain one grid.
	//CENTER: we assume that the nodes are at cell centers
	//CORNER: we assume that the nodes are at cell corners
	enum GridType { CENTER = 0, CORNER };

	template<int d>
	class GridIndexer {
		Typedef_VectorD(d);
	public:
		//number of blocks along x,y,z axis
		int nbx, nby, nbz;
		//for d==2: bits 012, 345 are the padded length of dimension x,y
		//for d==3: bits 01, 23, 45 are the padded length of dimension x,y,z
		unsigned char padding;
		GridIndexer() :nbx(0), nby(0), nbz(0), padding(0) {}
		GridIndexer(const VectorDi _counts) {
			constexpr int block_size = Block_Size();
			VectorDi padded_counts = MathFunc::Round_Up_To_Align<d>(_counts, block_size);
			VectorDi diff_counts = padded_counts - _counts;
			if constexpr (d == 2) padding = diff_counts[0] | (diff_counts[1] << 3);
			else if constexpr (d == 3) padding = diff_counts[0] | (diff_counts[1] << 2) | (diff_counts[2] << 4);
			if constexpr (d == 2) {
				nbx = padded_counts[0] >> 3;
				nby = padded_counts[1] >> 3;
				nbz = 0;
			}
			else if constexpr (d == 3) {
				nbx = padded_counts[0] >> 2;
				nby = padded_counts[1] >> 2;
				nbz = padded_counts[2] >> 2;
			}
		}

		//============================part 1: nodes and general data===============================

		__host__ __device__ static constexpr int Block_Size(void) { return (d == 2 ? 8 : 4); }
		__host__ __device__ constexpr bool Is_Unpadded(void)const { return padding == 0; }
		__host__ __device__ constexpr VectorDi Memory_Counts(void) const {
			if constexpr (d == 2) return Vector2i(nbx << 3, nby << 3);
			else if constexpr (d == 3) return Vector3i(nbx << 2, nby << 2, nbz << 2);
			else Assert(false, "Grid::Index: d==2 or d==3");
		}
		__host__ __device__ constexpr int Memory_Size(void)const { return Memory_Counts().prod(); }
		__host__ __device__ constexpr int Dimension(const int axis) const {
			if constexpr (d == 2) {
				if (axis == 0) return (nbx << 3) - (padding & 7);
				else if (axis == 1) return (nby << 3) - ((padding >> 3) & 7);
			}
			else if constexpr (d == 3) {
				if (axis == 0) return (nbx << 2) - (padding & 3);
				else if (axis == 1) return (nby << 2) - ((padding >> 2) & 3);
				else if (axis == 2) return (nbz << 2) - ((padding >> 4) & 3);
			}
		}
		__host__ __device__ constexpr VectorDi Counts(void)const {
			if constexpr (d == 2) return Vector2i(Dimension(0), Dimension(1));
			else if constexpr (d == 3) return Vector3i(Dimension(0), Dimension(1), Dimension(2));
		}
		__host__ __device__ constexpr int Valid_Size(void)const { return Counts().prod(); }
		__host__ __device__ constexpr bool Valid(const int i, const int j = 0, const int k = 0)const {
			if constexpr (d == 2) return 0 <= i && i < Dimension(0) && 0 <= j && j < Dimension(1);
			else if constexpr (d == 3) return 0 <= i && i < Dimension(0) && 0 <= j && j < Dimension(1) && 0 <= k && k < Dimension(2);
		}
		__host__ __device__ constexpr bool Valid(const VectorDi coord)const {
			if constexpr (d == 2) return Valid(coord[0], coord[1]);
			else if constexpr (d == 3) return Valid(coord[0], coord[1], coord[2]);
		}
		//__host__ __device__ constexpr bool Valid_Index(const int index)const { return Valid(Coord(index)); }

		__host__ __device__ constexpr int Index(const int i, const int j = 0, const int k = 0) const {
			if constexpr (d == 2) {
				//y-x in each block
				return ((j >> 3) * nbx + (i >> 3)) * 64 + ((j & 7) * 8 + (i & 7));
			}
			else if constexpr (d == 3) {
				int bx = i >> 2, by = j >> 2, bz = k >> 2;
				int idx = i & 0b11, idy = j & 0b11, idz = k & 0b11;
				//z-y-x in each block
				return ((bz * nby + by) * nbx + bx) * 64 + ((idz * 4 + idy) * 4 + idx);
			}
			else Assert(false, "Grid::Index: d==2 or d==3");
		}

		__host__ __device__ constexpr int Index(const VectorDi coord) const {
			if constexpr (d == 2) return Index(coord[0], coord[1]);
			else if constexpr (d == 3) return Index(coord[0], coord[1], coord[2]);
			else Assert(false, "Grid::Index: d==2 or d==3");
		}

		__host__ __device__ constexpr VectorDi Coord(const int index)const {
			if constexpr (d == 2) {
				int idx = index & 0b111;
				int idy = (index & 0b111000) >> 3;
				int b = index >> 6;
				int x = ((b % nbx) << 3) + idx;
				int y = ((b / nbx) << 3) + idy;
				return Vector2i(x, y);
			}
			else if constexpr (d == 3) {
				int idx = index & 0b11, idy = (index & 0b1100) >> 2, idz = (index & 0b110000) >> 4;
				int b = index >> 6;
				int x = ((b % nbx) << 2) + idx;
				int y = (((b / nbx) % nby) << 2) + idy;
				int z = ((b / nbx / nby) << 2) + idz;
				return Vector3i(x, y, z);
			}
		}

		//============================part 2: syntactic sugar for faces===============================

		__host__ __device__ constexpr int Valid_Face(const int axis, const VectorDi face)const {
			VectorDi face_counts = Counts(); face_counts[axis]++;
			if constexpr (d == 2) return 0 <= face[0] && face[0] < face_counts[0] && 0 <= face[1] && face[1] < face_counts[1];
			else if constexpr (d == 3) return 0 <= face[0] && face[0] < face_counts[0] && 0 <= face[1] && face[1] < face_counts[1] && 0 <= face[2] && face[2] < face_counts[2];
		}
		__host__ __device__ constexpr int Face_Index(const int axis, const VectorDi face) const
		{
			if constexpr (d == 2) {
				int x = face[0], y = face[1];
				int b_ind = (y >> 3) * (nbx + (axis == 0 && !(padding & 7))) + (x >> 3);
				int idx = x & 7, idy = y & 7;
				return b_ind * 64 + (idy * 8 + idx);
			}
			else if constexpr (d == 3) {
				int x = face[0], y = face[1], z = face[2];
				int fnbx = nbx + (axis == 0 && !(padding & 3)), fnby = nby + (axis == 1 && !((padding >> 2) & 3));
				int bx = x >> 2, by = y >> 2, bz = z >> 2;
				int idx = x & 0b11, idy = y & 0b11, idz = z & 0b11;
				return ((bz * fnby + by) * fnbx + bx) * 64 + ((idz * 4 + idy) * 4 + idx);
			}
		}
		__host__ __device__ constexpr VectorDi Face_Coord(const int axis, const int face_index)const {
			VectorDi face;
			if constexpr (d == 2) {
				int b_ind = face_index >> 6;
				int i1 = (face_index & 0b111000) >> 3, i0 = face_index & 0b111;
				int fnbx = nbx + (axis == 0);
				int idx, idy;

				idx = i0, idy = i1;

				face[0] = ((b_ind % fnbx) << 3) + idx;
				face[1] = ((b_ind / fnbx) << 3) + idy;
			}
			else if constexpr (d == 3) {
				int fnbx = nbx + (axis == 0), fnby = nby + (axis == 1);
				int idx = face_index & 0b11, idy = (face_index & 0b1100) >> 2, idz = (face_index & 0b110000) >> 4;
				int b = face_index >> 6;

				face[0] = ((b % fnbx) << 2) + idx;
				face[1] = (((b / fnbx) % fnby) << 2) + idy;
				face[2] = ((b / fnbx / fnby) << 2) + idz;
			}
			return face;
		}

		//============================part 3: iterators===============================

		std::tuple<dim3, dim3> Get_Kernel_Dims(void) const {
			dim3 blocknum, blocksize;
			if constexpr (d == 2) {
				blocknum = dim3(nbx, nby);
				blocksize = dim3(8, 8);
			}
			else if constexpr (d == 3) {
				blocknum = dim3(nbx, nby, nbz);
				blocksize = dim3(4, 4, 4);
			}
			return std::make_tuple(blocknum, blocksize);
		}
		//NOTE: WILL EXECUTE ALL THE MEMORY COUNTS INSTEAD OF ONLY VALID CELLS
		//MUST CHECK IF VALID INSIDE THE KERNEL FUNCTION
		template<class Func, class ...Args>
		void Exec_Kernel(Func kernel_func, const Args&...args) const {
			auto [blocknum, blocksize] = Get_Kernel_Dims();
			kernel_func << <blocknum, blocksize >> > (args...);
		}
	};

	template<int d>
	constexpr bool operator == (const GridIndexer<d>&& grid1, decltype(grid1) grid2) {
		return grid1.nbx == grid2.nbx && grid1.nby == grid2.nby && grid1.nbz == grid2.nbz && grid1.padding == grid2.padding;
	}

	template<int d>
	class Grid: public GridIndexer<d> {
		Typedef_VectorD(d);
		using Base = GridIndexer<d>;
	public:
		VectorD pos_min;
		real dx;
		

		Grid(const VectorDi _counts = VectorDi::Zero(), const real _dx = 0, const VectorD domain_min = VectorD::Zero(), const GridType gtype = CENTER) :
			GridIndexer<d>(_counts),
			dx(_dx)
		{
			if (gtype == CORNER) pos_min = domain_min;
			else pos_min = domain_min + MathFunc::V<d>(0.5, 0.5, 0.5) * dx;
		}

		constexpr GridIndexer<d> Indexer(void)const { return (*this); }

//=============================================First part: basic data==========================================================
		__host__ __device__ real Dx(void) { return dx; }

		__host__ __device__ constexpr VectorD Position(const VectorDi node)const {
			return pos_min + node.template cast<real>() * dx;
		}
		__host__ __device__ constexpr VectorD Node_Min(void)const { return pos_min; }
		__host__ __device__ constexpr VectorD Node_Max(void)const { return Position(Base::Counts() - VectorDi::Ones()); }
		__host__ __device__ VectorD Center(void) const {
			//grid type doesn't matter
			return (real)0.5 * (Node_Min() + Node_Max());
		}
		__host__ __device__ VectorD Domain_Min(const GridType gtype) const {
			if (gtype == CENTER) {
				return pos_min - VectorD::Ones() * 0.5 * dx;
			}
			return pos_min;
		}
		__host__ __device__ VectorD Domain_Max(const GridType gtype) const {
			VectorD node_max = Node_Max();
			if (gtype == CENTER) {
				return node_max + VectorD::Ones() * 0.5 * dx;
			}
			return node_max;
		}

		__host__ __device__ void Get_Fraction(const VectorD pos, VectorDi& node, VectorD& frac)const {
			VectorD coord_with_frac = (pos - pos_min) / dx;
			node = coord_with_frac.template cast<int>();
			for (int i = 0; i < d; i++)frac[i] = coord_with_frac[i] - node[i];
		}



//========================================Second part: geometry operation=================================================================
		//// Navigate on adjacent nodes
		//added by Zhiqi Li
		__host__ __device__ constexpr static int Neighbor_Node_Number(void) { return d * 2; }
		__host__ __device__ static VectorDi Neighbor_Node(const VectorDi coord, const int i) {
			if constexpr (d == 2) {
				const Vector2i neighbor_offset[4] = { Vector2i(-1,0),Vector2i(0,-1),Vector2i(0,1),Vector2i(1,0) };
				return coord + neighbor_offset[i];
			}
			else if constexpr (d == 3) {
				const Vector3i neighbor_offset[6] = { Vector3i(-1,0,0),Vector3i(0,-1,0),Vector3i(0,0,-1),Vector3i(0,0,1),Vector3i(0,1,0),Vector3i(1,0,0) };
				return coord + neighbor_offset[i];
			}
		}
		__host__ __device__ static VectorDi Neighbor_Node(const VectorDi coord, const int axis, const int side) {
			if constexpr (d == 2) {
				const int nb_ids[2][2] = { 0,3,1,2 }; return Grid<d>::Neighbor_Node(coord, nb_ids[axis][side]);
			}
			else if constexpr (d == 3) {
				const int nb_ids[3][2] = { 0,5,1,4,2,3 }; return Grid<d>::Neighbor_Node(coord, nb_ids[axis][side]);
			}
		}
		__host__ __device__ static int Neighbor_Node_Axis(const int i) {
			if constexpr (d == 2) {
				const int axies[4] = { 0,1,1,0 }; return axies[i];
			}
			else if constexpr (d == 3) {
				const int axies[6] = { 0,1,2,2,1,0 }; return axies[i];
			}
		}

		//Navigate on neighbor ring
		__host__ __device__ constexpr static int Neighbor_Ring_Number(void) { return std::pow(3, d);}
		__host__ __device__ static VectorDi Neighbor_Ring_Node(const VectorDi coord, const int index) {
			if constexpr (d == 2) {
				assert(index >= 0 && index < 9); int i = index / 3; int j = index % 3; return coord + Vector2i(-1 + i, -1 + j);
			}
			else if constexpr (d == 3) {
				assert(index >= 0 && index < 27); int i = index / 9; int m = index % 9; int j = m / 3; int k = m % 3; return coord + Vector3i(-1 + i, -1 + j, -1 + k);
			}
		}

		////Geometry interfaces
		//Explain the grid as a corner grid, extract the center grid
		__host__ __device__ Grid<d> Cell_Grid(void)const {
			return Grid<d>(Base::Counts() - VectorDi::Ones(), dx, pos_min, CENTER);
		}

		////Explain the grid as a MAC grid, extract faces
		//__host__ __device__ VectorDi Face_Memory_Counts(const int axis)const { VectorDi fcounts = indexer.Counts(); fcounts[axis] += Block_Size(); return fcounts; }
		//__host__ __device__ int Face_Memory_Size(int axis)const { return Face_Memory_Counts(axis).prod(); }
		__host__ __device__ VectorD Face_Min(const int axis)const {
			VectorD offset = -VectorD::Unit(axis) * 0.5 * dx;
			return pos_min + offset;
		}
		__host__ __device__ constexpr VectorD Face_Center(const int axis, const VectorDi face) const {
			return Face_Min(axis) + face.template cast<real>() * dx;
		}
		__host__ __device__ Grid<d> Face_Grid(const int axis)const {
			if (axis >= d) return Grid<d>();//placeholder
			VectorDi face_counts = Base::Counts(); face_counts[axis]++;
			return Grid<d>(face_counts, dx, Face_Min(axis), CORNER);
		}
//======================================================Third part: iterators==================================================================
		//serial iterator
		template<class CFunc>
		void Iterate_Nodes(CFunc f)const {
			const int memory_size = Base::Memory_Size();
			for (int c = 0; c < memory_size; c++) {
				const VectorDi node = Base::Coord(c);
				if (Base::Valid(node)) f(node);
			}
		}

		//parallel iterator
		template<class CFunc>
		void Exec_Nodes(CFunc f) const {
			const int memory_size = Base::Memory_Size();
#pragma omp parallel for
			for (int c = 0; c < memory_size; c++) {
				const VectorDi node = Base::Coord(c);
				if (Base::Valid(node))f(node);
			}
		}
		
		template<class ICFunc>
		void Iterate_Faces(ICFunc f) {
			for (int axis = 0; axis < d; axis++) {
				Grid<d> face_grid = Face_Grid(axis);
				auto face_node_f = std::bind(f, axis, std::placeholders::_1);
				face_grid.Iterate_Nodes(face_node_f);
			}
		}

		template<class ICFunc>
		void Exec_Faces(ICFunc f) const {
			for (int axis = 0; axis < d; axis++) {
				Grid<d> face_grid = Face_Grid(axis);
				auto face_node_f = std::bind(f, axis, std::placeholders::_1);
				face_grid.Exec_Nodes(face_node_f);
			}
		}
	};

}


//fmt adaptor for GridIndexer
template <int d>
struct fmt::formatter<Meso::GridIndexer<d>> {
	constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin()) {
		//https://fmt.dev/latest/api.html#udt
		auto it = ctx.begin(), end = ctx.end();
		if (it != end && *it != '}') throw format_error("invalid format");
		// Return an iterator past the end of the parsed range:
		return it;
	}

	// Formats the point p using the parsed format specification (presentation)
	// stored in this formatter.
	template <typename FormatContext>
	auto format(const Meso::GridIndexer<d> G, FormatContext& ctx) -> decltype(ctx.out()) {
		std::string out = "";
		Meso::Vector<int, d> valid_counts = G.Counts();
		Meso::Vector<int, d> memory_counts = G.Memory_Counts();
		if constexpr (d == 2) out += fmt::format("[{},{}]->[{},{}]", valid_counts[0], valid_counts[1], memory_counts[0], memory_counts[1]);
		else if constexpr (d == 3)out += fmt::format("[{},{},{}]->[{},{},{}]", valid_counts[0], valid_counts[1], valid_counts[2], memory_counts[0], memory_counts[1], memory_counts[2]);
		return format_to(ctx.out(), "{}", out);
	}
};

//fmt adaptor for GridIndexer
template <int d>
struct fmt::formatter<Meso::Grid<d>> {
	constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin()) {
		//https://fmt.dev/latest/api.html#udt
		auto it = ctx.begin(), end = ctx.end();
		if (it != end && *it != '}') throw format_error("invalid format");
		// Return an iterator past the end of the parsed range:
		return it;
	}

	// Formats the point p using the parsed format specification (presentation)
	// stored in this formatter.
	template <typename FormatContext>
	auto format(const Meso::Grid<d> G, FormatContext& ctx) -> decltype(ctx.out()) {
		std::string out = "";
		Meso::Vector<int, d> valid_counts = G.Counts();
		Meso::Vector<int, d> memory_counts = G.Memory_Counts();
		if constexpr (d == 2) out += fmt::format("[{},{}]->[{},{}] node box {}-{}", valid_counts[0], valid_counts[1], memory_counts[0], memory_counts[1], G.Node_Min(), G.Node_Max());
		else if constexpr (d == 3)out += fmt::format("[{},{},{}]->[{},{},{}] node box {}-{}", valid_counts[0], valid_counts[1], valid_counts[2], memory_counts[0], memory_counts[1], memory_counts[2], G.Node_Min(), G.Node_Max());
		return format_to(ctx.out(), "{}", out);
	}
};
