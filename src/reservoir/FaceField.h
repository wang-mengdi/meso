//////////////////////////////////////////////////////////////////////////
// Face data on MacGrid
// Copyright (c) (2022-), Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "Grid.h"
#include "Field.h"

namespace Meso {

	template<class T, int d, DataHolder side = DataHolder::HOST>
	class FaceField {
		Typedef_VectorD(d);
	public:
		Grid<d> grid;
		//std::array<Array<T, side>, d> face_data;
		std::shared_ptr<Array<T, side>> face_data[3] = { nullptr,nullptr,nullptr };
		FaceField() {}
		FaceField(const Grid<d>& _grid) { Init(_grid); }
		FaceField(const Grid<d>& _grid, const T value) { Init(_grid);  Fill(value); }
		template<DataHolder side1> FaceField(const FaceField<T, d, side1>& f1) { Deep_Copy(f1); }
		template<DataHolder side1> FaceField(FaceField<T, d, side1>& f1) { Deep_Copy(f1); }
		void Fill(const T value) { for (int axis = 0; axis < d; axis++) ArrayFunc::Fill(*face_data[axis], value); }
		void Init(const Grid<d>& _grid, const T value) { Init(_grid); Fill(value); }
		void Init(const Grid<d>& _grid) {
			grid = _grid;
			for (int axis = 0; axis < d; axis++) {
				Grid<d> face_grid = grid.Face_Grid(axis);
				int n = face_grid.Memory_Size();
				if (face_data[axis] == nullptr) face_data[axis] = std::make_shared<Array<T, side>>(n);
				else face_data[axis]->resize(n);
				checkCudaErrors(cudaGetLastError());
			}
		}
		template<DataHolder side1>
		void Init(const FaceField<T, d, side1>& f1) {
			Deep_Copy(f1);
		}
		bool Empty(void)const { for (int axis = 0; axis < d; axis++) if (face_data[axis] == 0) return true; return false; }
		template<DataHolder side1> 
		void Deep_Copy(const FaceField<T, d, side1>& f1) {
			Init(f1.grid);
			for (int i = 0; i < d; i++) {
				//deep copy
				ArrayFunc::Copy(*face_data[i], f1.Data(i));
			}
		}

		template<DataHolder side1>
		FaceField<T, d, side>& operator = (const FaceField<T, d, side1>& f1) { Deep_Copy(f1); return *this; }
		template<DataHolder side1>
		FaceField<T, d, side>& operator = (FaceField<T, d, side1>& f1) { Deep_Copy(f1); return *this; }
		inline T& operator()(const int axis, const VectorDi face) { return (*(face_data[axis]))[grid.Face_Index(axis, face)]; }
		inline const T& operator()(int axis, const VectorDi face) const { return (*(face_data[axis]))[grid.Face_Index(axis, face)]; }
		inline const T Get(int axis, const VectorDi face) const { return (*(face_data[axis]))[grid.Face_Index(axis, face)]; }

		constexpr Array<T, side>& Data(const int axis)noexcept { return *face_data[axis]; }
		constexpr const Array<T, side>& Data(const int axis)const noexcept { return *face_data[axis]; }
		constexpr T* Data_Ptr(const int axis) noexcept { return face_data[axis] == nullptr ? nullptr : thrust::raw_pointer_cast(face_data[axis]->data()); }
		constexpr const T* Data_Ptr(const int axis) const noexcept {
			return face_data[axis] == nullptr ? nullptr : thrust::raw_pointer_cast(face_data[axis]->data());
		}
		constexpr Field<T, d, side> Face_Reference(const int axis)const {
			Assert(face_data[axis] != nullptr, "Field::Face_Reference error: empty face_data[{}]", axis);
			//will reference original data
			return Field<T, d, side>(grid.Face_Grid(axis), face_data[axis]);
		}

		void operator += (const Vector<T, d> vec) {
			for (int axis = 0; axis < d; axis++) {
				ArrayFunc::Add_Scalar(Data(axis), vec[axis]);
			}
		}
		void operator += (const FaceField<T, d, side>& f1) {
			for (int axis = 0; axis < d; axis++) {
				ArrayFunc::Add(Data(axis), f1.Data(axis));
			}
		}
		void operator -= (const FaceField<T, d, side>& f1) {
			for (int axis = 0; axis < d; axis++) {
				ArrayFunc::Minus(Data(axis), f1.Data(axis));
			}
		}
		void operator *= (const FaceField<T, d, side>& f1) {
			for (int axis = 0; axis < d; axis++) {
				ArrayFunc::Multiply(Data(axis), f1.Data(axis));
			}
		}
		template<class T1> void operator *= (const T1 a) {
			for (int axis = 0; axis < d; axis++) ArrayFunc::Multiply_Scalar(Data(axis), a);
		}
	
		template<class ICFunc>
		void Iterate_Faces(ICFunc f) {
			grid.Iterate_Faces(f);
		}

		template<class ICFunc>
		void Exec_Faces(ICFunc f) {
			grid.Exec_Faces(f);
		}

		template<class ICFuncT>
		void Calc_Faces(ICFuncT f) {
			Grid<d> grid2 = grid;
			for (int axis = 0; axis < d; axis++) {
				auto face_node_f = std::bind(f, axis, std::placeholders::_1);
				Field<T, d, side> face_field = Face_Reference(axis);
				face_field.Calc_Nodes(face_node_f);
			}
		}
	};


	template<class T, int d> using FaceFieldDv = FaceField<T, d, DEVICE>;

}

//fmt adaptor for FaceField
template <class T, Meso::DataHolder side>
struct fmt::formatter<Meso::FaceField<T, 2, side>> {
	constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin()) {
		//https://fmt.dev/latest/api.html#udt
		auto it = ctx.begin(), end = ctx.end();
		if (it != end && *it != '}') throw format_error("invalid format");
		// Return an iterator past the end of the parsed range:
		return it;
	}

	void Append_String(const Meso::FaceField<T, 2>& F, std::string& out) {
		fmt::formatter<Meso::Field<T, 2, side>> fmt1;
		for (int axis = 0; axis < 2; axis++) {
			out += "axis: "; out += std::to_string(axis); out += "\n";
			Meso::Field<T, 2> fi = F.Face_Reference(axis);
			fmt1.Append_String(fi, out);
			out += "===========\n";
		}
	}

	// Formats the point p using the parsed format specification (presentation)
	// stored in this formatter.
	template <typename FormatContext>
	auto format(const Meso::FaceField<T, 2, side>& F, FormatContext& ctx) -> decltype(ctx.out()) {
		std::string out;
		if constexpr (side == Meso::DataHolder::HOST) Append_String(F, out);
		else if constexpr (side == Meso::DataHolder::DEVICE) {
			Meso::FaceField<T, 2> F_host = F;
			Append_String(F_host, out);
		}
		return format_to(ctx.out(), "{}", out);
	}
};

template<class T, Meso::DataHolder side>
struct fmt::formatter<Meso::FaceField<T, 3, side>> {
	constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin()) {
		//https://fmt.dev/latest/api.html#udt
		auto it = ctx.begin(), end = ctx.end();
		if (it != end && *it != '}') throw format_error("invalid format");
		// Return an iterator past the end of the parsed range:
		return it;
	}

	void Update_String(const Meso::FaceField<T, 3>& F, std::string& out) {
		fmt::formatter<Meso::Field<T, 3, side>> fmt1;
		for (int axis = 0; axis < 3; axis++) {
			out += "axis: "; out += std::to_string(axis); out += "\n";
			Meso::Field<T, 3> fi = F.Face_Reference(axis);
			fmt1.Append_String(fi, out);
			out += "===========\n";
		}
	}

	// Formats the point p using the parsed format specification (presentation)
	// stored in this formatter.
	template <typename FormatContext>
	auto format(const Meso::FaceField<T, 3, side>& F, FormatContext& ctx) -> decltype(ctx.out()) {
		std::string out = "";
		if constexpr (side == Meso::DataHolder::HOST) Update_String(F, out);
		else {
			Meso::FaceField<T, 3> F_host = F;
			Update_String(F_host, out);
		}
		return format_to(ctx.out(), "{}", out);
	}
};
