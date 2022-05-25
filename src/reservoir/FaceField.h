//////////////////////////////////////////////////////////////////////////
// Face data on MacGrid
// Copyright (c) (2022-), Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "Grid.h"


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
		template<DataHolder side1> FaceField(const FaceField<T, d, side1>& f1) { Init(f1); }
		void Fill(const T value) { for (int axis = 0; axis < d; axis++) ArrayFunc::Fill(*face_data[axis], value); }
		void Init(const Grid<d>& _grid, const T value) { Init(_grid); Fill(value); }
		void Init(const Grid<d>& _grid) {
			grid = _grid;
			for (int axis = 0; axis < d; axis++) {
				int n = grid.Face_DoF(axis);
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
				*face_data[i] = f1.Data(i);
			}
		}

		template<DataHolder side1>
		FaceField<T, d, side>& operator = (const FaceField<T, d, side1>& f1) {
			Deep_Copy(f1);
			return *this;
		}
		template<DataHolder side1>
		FaceField<T, d, side>& operator = (FaceField<T, d, side1>& f1) {
			Deep_Copy(f1);
			return *this;
		}
		inline T& operator()(const int axis, const VectorDi face) { return (*(face_data[axis]))[grid.Face_Index(axis, face)]; }
		inline const T& operator()(int axis, const VectorDi face) const { return (*(face_data[axis]))[grid.Face_Index(axis, face)]; }
		inline const T Get(int axis, const VectorDi face) const { return (*(face_data[axis]))[grid.Face_Index(axis, face)]; }

		constexpr Array<T, side>& Data(const int axis)noexcept { return *face_data[axis]; }
		constexpr const Array<T, side>& Data(const int axis)const noexcept { return *face_data[axis]; }
		constexpr T* Data_Ptr(const int axis) noexcept { return face_data[axis] == nullptr ? nullptr : thrust::raw_pointer_cast(face_data[axis]->data()); }
		constexpr const T* Data_Ptr(const int axis) const noexcept {
			return face_data[axis] == nullptr ? nullptr : thrust::raw_pointer_cast(face_data[axis]->data());
		}

		void operator += (const Vector<T, d> vec) {
			for (int axis = 0; axis < d; axis++) {
				ArrayFunc::Add(Data(axis), vec[axis]);
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

		T Max_Abs(void) {
			real max_val = 0;
			for (int axis = 0; axis < d; axis++) {
				max_val = std::max<T>(max_val, ArrayFunc::Max_Abs<T>(Data(axis)));
			}
			return max_val;
		}

		template<class IFFunc>
		void Iterate_Faces(IFFunc f) {
			for (int axis = 0; axis < d; axis++) {
				int n = grid.Face_DoF(axis);
				for (int i = 0; i < n; i++) {
					VectorDi face = grid.Face_Coord(axis, i);
					f(axis, face);
				}
			}
		}

		template<class IFFuncT>
		void Calc_Faces(IFFuncT f) {
			Grid<d> grid2 = grid;
			for (int axis = 0; axis < d; axis++) {
				Assert(face_data[axis] != nullptr, "FaceField::Calc_Faces error: nullptr data at axis {}", axis);
				const int dof = grid.Face_DoF(axis);
				/// Note!
				/// Why I write the lambda out here:
				/// A strange error: For this host platform, an extended lambda cannot be defined inside the 'if'
				/// or 'else' block of a constexpr if statement
				if constexpr (side == DEVICE) {
					auto f_device = [f, axis, grid2]__device__(const int idx) {
						return f(axis, grid2.Face_Coord(axis, idx));
					};
					thrust::counting_iterator<int> idxfirst(0);
					thrust::counting_iterator<int> idxlast = idxfirst + dof;
					thrust::transform(
						idxfirst,
						idxlast,
						face_data[axis]->begin(),
						f_device
					);
				}
				else {
#pragma omp parallel for
					for (int i = 0; i < dof; i++) {
						VectorDi face = grid2.Face_Coord(axis, i);
						(*this)(axis, face) = f(axis, face);
					}
				}
			}
			if constexpr (side == DEVICE) {
				cudaFree(grid_gpu);
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

	void Update_String(const Meso::FaceField<T, 2>& F, std::string& out) {
		out = "";
		for (int axis = 0; axis < 2; axis++) {
			out += "axis: "; out += std::to_string(axis); out += "\n";
			for (int i = 0; i < F.grid.counts[0]; i++) {
				for (int j = 0; j < F.grid.counts[1]; j++) {
					out += Meso::StringFunc::To_String_Simple(F(axis,Eigen::Vector2i(i, j))) + " ";
				}
				out += "\n";
			}
			out += "===========\n";
		}
	}

	// Formats the point p using the parsed format specification (presentation)
	// stored in this formatter.
	template <typename FormatContext>
	auto format(const Meso::FaceField<T, 2, side>& F, FormatContext& ctx) -> decltype(ctx.out()) {
		std::string out;
		if constexpr (side == Meso::DataHolder::HOST) Update_String(F, out);
		else if constexpr (side == Meso::DataHolder::DEVICE) {
			Meso::FaceField<T, 2> F_host = F;
			Update_String(F_host, out);
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
		out = "";
		for (int axis = 0; axis < 3; axis++) {
			out += "axis: "; out += std::to_string(axis); out += "\n";
			for (int i = 0; i < F.grid.counts[0]; i++) {
				for (int j = 0; j < F.grid.counts[1]; j++) {
					for (int k = 0; k < F.grid.counts[2]; k++) {
						out += Meso::StringFunc::To_String_Simple(F(axis,Eigen::Vector3i(i, j, k))) + " ";
					}
					out += "\n";
				}
				out += "===========\n";
			}
			out += "======================\n";
		}
	}

	// Formats the point p using the parsed format specification (presentation)
	// stored in this formatter.
	template <typename FormatContext>
	auto format(const Meso::FaceField<T, 3, side>& F, FormatContext& ctx) -> decltype(ctx.out()) {
		std::string out;
		if constexpr (side == Meso::DataHolder::HOST) Update_String(F, out);
		else {
			Meso::FaceField<T, 3> F_host = F;
			Update_String(F_host, out);
		}
		return format_to(ctx.out(), "{}", out);
	}
};