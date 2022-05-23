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

	template<class T, int d, DataHolder side = HOST>
	class Field {
		Typedef_VectorD(d);
	public:
		Grid<d> grid;
		std::shared_ptr<Array<T, side>> data = nullptr;
		Field() {}
		Field(const Grid<d>& _grid, std::shared_ptr<Array<T, side>> _data) { grid = _grid; data = _data; }
		Field(const Grid<d>& _grid) { Init(_grid); }
		Field(const Grid<d> _grid, const T value) { Init(_grid, value); }
		template<DataHolder side1> Field(const Field<T, d, side1>& f1) { *this = f1; }
		void Init(const Grid<d> _grid) {
			grid = _grid;
			if (data == nullptr) data = std::make_shared<Array<T, side>>(grid.DoF());
			else data->resize(grid.DoF());
			//data.resize(grid.DoF());
			checkCudaErrors(cudaGetLastError());
		}
		void Init(const Grid<d> _grid, const T value) {
			Init(_grid);
			ArrayFunc::Fill(*data, value);
		}
		void Fill(const T value) { ArrayFunc::Fill(Data(), value); }
		bool Empty(void)const { return data == nullptr; }
		constexpr Array<T, side>& Data(void) noexcept {
			return *data;
		}
		constexpr const Array<T, side>& Data(void)const noexcept {
			return *data;
		}
		constexpr T* Data_Ptr(void) noexcept {
			return thrust::raw_pointer_cast(data->data());
		}
		constexpr const T* Data_Ptr(void) const noexcept {
			return thrust::raw_pointer_cast(data->data());
		}

		template<DataHolder side1>
		void Deep_Copy(const Field<T, d, side1>& f1) {
			Init(f1.grid);
			//deep copy
			*data = f1.Data();
		}


		template<DataHolder side1>
		Field<T, d, side>& operator = (const Field<T, d, side1>& f1) {
			Deep_Copy(f1);
			return *this;
		}
		template<DataHolder side1>
		Field<T, d, side>& operator = (Field<T, d, side1>& f1) {
			Deep_Copy(f1);
			return *this;
		}
		inline T& operator()(const VectorDi coord) { return (*data)[grid.Index(coord)]; }
		inline const T& operator()(const VectorDi coord) const { return (*data)[grid.Index(coord)]; }
		const T Get(const VectorDi coord)const { return (*data)[grid.Index(coord)]; }
		template<class T1> void operator *= (const T1 a) { ArrayFunc::Multiply_Scalar(Data(), a); }
		void operator -= (const Field<T, d, side>& f1) { ArrayFunc::Minus(Data(), f1.Data()); }

		T Max_Abs(void) {
			return ArrayFunc::Max_Abs<T>(Data());
		}

		template<class CFunc>
		void Iterate_Nodes(CFunc f) {
			const int dof = grid.DoF();
			for (int c = 0; c < dof; c++) {
				f(grid.Coord(c));
			}
		}

		/// Modify by Zhiqi Li, add surpoort for GPU
		template<class CFuncT>//CFuncT is a function: VectorDi->T, takes the cell index
		__host__ void Calc_Nodes(CFuncT f) {
			Assert(data != nullptr, "Field::Calc_Cells error: nullptr data");
			const int dof = grid.DoF();
			thrust::counting_iterator<int> idxfirst(0);
			thrust::counting_iterator<int> idxlast = idxfirst + dof;
			Grid<d> grid2 = grid;
			thrust::transform(
				idxfirst,
				idxlast,
				data->begin(),
				[f, grid2]__device__ __host__(const int idx) {
				return f(grid2.Coord(idx));
			}
			);
		}

		template<class Fnode>//Fnode is a (void) function takes a node index
		void Exec_Nodes(Fnode f) const {
			return grid.Exec_Nodes(f);
		}

		template<class F, class ...Args>
		void Exec_Kernel(F kernel_func, const Args&...args) const {
			grid.Exec_Kernel(kernel_func, args...);
		}
	};

	template<class T, int d> using FieldDv = Field<T, d, DEVICE>;

}

//fmt adaptor for Field
template <class T, Meso::DataHolder side>
struct fmt::formatter<Meso::Field<T, 2, side>> {
	constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin()) {
		//https://fmt.dev/latest/api.html#udt
		auto it = ctx.begin(), end = ctx.end();
		if (it != end && *it != '}') throw format_error("invalid format");
		// Return an iterator past the end of the parsed range:
		return it;
	}

	void Update_String(const Meso::Field<T, 2>& F, std::string& out) {
		out = "";
		//out += to_string(F.grid.counts[0]);
		for (int i = 0; i < F.grid.counts[0]; i++) {
			for (int j = 0; j < F.grid.counts[1]; j++) {
				out += Meso::StringFunc::To_String_Simple(F(Eigen::Vector2i(i, j))) + " ";
			}
			out += "\n";
		}
	}

	// Formats the point p using the parsed format specification (presentation)
	// stored in this formatter.
	template <typename FormatContext>
	auto format(const Meso::Field<T, 2, side>& F, FormatContext& ctx) -> decltype(ctx.out()) {
		std::string out;
		if constexpr (side == Meso::DataHolder::HOST) Update_String(F, out);
		else if constexpr (side == Meso::DataHolder::DEVICE) {
			Meso::Field<T, 2> F_host = F;
			Update_String(F_host, out);
		}
		return format_to(ctx.out(), "{}", out);
	}
};

template<class T, Meso::DataHolder side>
struct fmt::formatter<Meso::Field<T, 3, side>> {
	constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin()) {
		//https://fmt.dev/latest/api.html#udt
		auto it = ctx.begin(), end = ctx.end();
		if (it != end && *it != '}') throw format_error("invalid format");
		// Return an iterator past the end of the parsed range:
		return it;
	}

	void Update_String(const Meso::Field<T, 3>& F, std::string& out) {
		out = "";
		//out += to_string(F.grid.counts[0]);
		for (int i = 0; i < F.grid.counts[0]; i++) {
			for (int j = 0; j < F.grid.counts[1]; j++) {
				for (int k = 0; k < F.grid.counts[2]; k++) {
					out += Meso::StringFunc::To_String_Simple(F(Eigen::Vector3i(i, j, k))) + " ";
				}
				out += "\n";
			}
			out += "===========\n";
		}
	}

	// Formats the point p using the parsed format specification (presentation)
	// stored in this formatter.
	template <typename FormatContext>
	auto format(const Meso::Field<T, 3, side>& F, FormatContext& ctx) -> decltype(ctx.out()) {
		std::string out;
		if constexpr (side == Meso::DataHolder::HOST) Update_String(F, out);
		else {
			Meso::Field<T, 3> F_host = F;
			Update_String(F_host, out);
		}
		return format_to(ctx.out(), "{}", out);
	}
};