//////////////////////////////////////////////////////////////////////////
// Basic grid data (with data included)
// Copyright (c) (2022-), Bo Zhu, Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "Grid.h"
#include <thrust/iterator/zip_iterator.h>
#include <thrust/iterator/counting_iterator.h>

namespace Meso {

	enum CellType { Air = 0, Fluid, Solid };

	template<class T, int d, DataHolder side = DataHolder::HOST>
	class Field {
		Typedef_VectorD(d);
	public:
		Grid<d, GridType::CENTER> grid;
		Array<T, side> data;
		Field() {}
		Field(const Grid<d, GridType::CENTER>& _grid) { Init(_grid); }
		Field(const Grid<d, GridType::CENTER> _grid, const T value) { Init(_grid, value); }
		template<DataHolder side1> Field(const Field<T, d, side1>& f1) { *this = f1; }
		void Init(const Grid<d, GridType::CENTER> _grid) {
			grid = _grid;
			data.resize(grid.DoF());
		}
		void Init(const Grid<d, GridType::CENTER> _grid, const T value) {
			Init(_grid);
			ArrayFunc::Fill(data, value);
		}
		template<DataHolder side1>
		Field<T, d, side>& operator = (const Field<T, d, side1>& f1) {
			grid = f1.grid;
			data = f1.data;
			return *this;
		}

		constexpr T* Data(void) noexcept {
			if constexpr (side == HOST) return data.data();
			else return thrust::raw_pointer_cast(data.data());
		}
		constexpr const T* Data(void) const noexcept {
			if constexpr (side == HOST) return data.data();
			else return thrust::raw_pointer_cast(data.data());
		}
		
		template<DataHolder side1> void Copy(const Field<T, d, side1>& f1) { Init(f1.grid); ArrayFunc::Copy(data, f1.data); }

		inline T& operator()(const VectorDi coord) { return data[grid.Index(coord)]; }
		inline const T& operator()(const VectorDi coord) const { return data[grid.Index(coord)]; }

		template<class CFunc>
		void Iterate_Cells(CFunc f){
			const int dof = grid.DoF();
			for (int c = 0; c < dof; c++) {
				f(grid.Coord(c));
			}
		}

		template<class CFuncT>//CFuncT is a function: VectorDi->T, takes the cell index
		void Calc_Cells(CFuncT f) {
			const int dof = grid.DoF();
			thrust::counting_iterator<int> idxfirst(0);
			thrust::counting_iterator<int> idxlast = idxfirst + dof;
			thrust::transform(
				idxfirst,
				idxlast,
				data.begin(),
				[f, this](const int idx) {
					return f(grid.Coord(idx));
				}
			);
		}
	};

	template<class T, int d> using FieldDv = Field<T, d, DEVICE>;



}

//fmt adaptor for Field
template <class T>
struct fmt::formatter<Meso::Field<T, 2>> {
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
	auto format(const Meso::Field<T, 2>& F, FormatContext& ctx) -> decltype(ctx.out()) {
		std::string out = "";
		//out += to_string(F.grid.counts[0]);
		for (int i = 0; i < F.grid.counts[0]; i++) {
			for (int j = 0; j < F.grid.counts[1]; j++) {
				out += Meso::IOFunc::To_String_Simple(F(Eigen::Vector2i(i, j))) + " ";
			}
			out += "\n";
		}
		return format_to(ctx.out(), "{}", out);
	}
};

template<class T>
struct fmt::formatter<Meso::Field<T, 3>> {
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
	auto format(const Meso::Field<T, 3>& F, FormatContext& ctx) -> decltype(ctx.out()) {
		std::string out = "";
		//out += to_string(F.grid.counts[0]);
		for (int i = 0; i < F.grid.counts[0]; i++) {
			for (int j = 0; j < F.grid.counts[1]; j++) {
				for (int k = 0; k < F.grid.counts[2]; k++) {
					out += Meso::IOFunc::To_String_Simple(F(Eigen::Vector3i(i, j, k))) + " ";
				}
				out += "\n";
			}
			out += "===========\n";
		}
		return format_to(ctx.out(), "{}", out);
	}
};