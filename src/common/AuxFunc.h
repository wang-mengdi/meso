//////////////////////////////////////////////////////////////////////////
// Common auxillary functions
// Copyright (c) (2022-), Bo Zhu, Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "Common.h"
#include <thrust/inner_product.h>
#include <thrust/transform.h>
#include <thrust/functional.h>
using namespace thrust::placeholders;

namespace Meso {

	namespace VectorFunc {
		////create vectors with compatible dimensions
		template<int d> Vector<real, d> V(const real x = (real)0, const real y = (real)0, const real z = (real)0);
		template<int d> Vector<int, d> Vi(const int x = 0, const int y = 0, const int z = 0, const int w = 0);
		template<int d> Vector<real, d> V(const Vector2 v2);
		template<int d> Vector<real, d> V(const Vector3 v2);

		////Round up vector to a multiple of bn
		template<int d> Vector<int, d> Round_Up_To_Align(Vector<int, d> v, int bn) {
			for (int i = 0; i < d; i++) {
				//round to the nearest multiple of block_size
				v[i] = ((v[i] + bn - 1) / bn) * bn;
			}
			return v;
		}
	}

	namespace ArrayFunc {
		template<class Array1, class T>
		void Fill(Array1& a, const T val) {
			thrust::fill(a.begin(), a.end(), val);
		}
		template<class T>
		T Dot(const Array<T, HOST>& a, decltype(a) b) {
			Assert(a.size() == b.size(), "[GPUFunc::Dot] try to dot length {} against {}", a.size(), b.size());
			return thrust::inner_product(a.begin(), a.end(), b.begin(), (T)0);
		}
		template<class T>
		T Dot(const ArrayDv<T>& a, decltype(a) b) {
			Assert(a.size() == b.size(), "[GPUFunc::Dot] try to dot length {} against {}", a.size(), b.size());
			return thrust::inner_product(a.begin(), a.end(), b.begin(), (T)0);
		}

		//element-wise multiplication, a*.=b
		template<class T, DataHolder side>
		void Multiply(Array<T, side>& a, const decltype(a) b) {
			thrust::transform(a.begin(), a.end(), b.begin(), a.begin(), thrust::multiplies<T>());
		}
		//element-wise add, a+.=b
		template<class T, DataHolder side>
		void Add(Array<T, side>& a, const decltype(a) b) {
			thrust::transform(a.begin(), a.end(), b.begin(), a.begin(), thrust::plus<T>());
		}
		//a=b, note it's reverse order of thrust::copy itself
		template<class Array1, class Array2>
		void Copy(Array1& a, const Array2& b) {
			thrust::copy(b.begin(), b.end(), a.begin());
		}
		//y=y+a*x
		template<class T>
		void Axpy(const real a, const ArrayDv<T>& x, ArrayDv<T>& y) {
			thrust::transform(x.begin(), x.end(), y.begin(), y.begin(), a * _1 + _2);
		}
		//x*=a
		template<class T>
		void Scal(const real a, ArrayDv<T>& x) {
			thrust::transform(x.begin(), x.end(), x.begin(), a * _1);
		}
		template<class TFuncT, class Array1, class Array2>
		void Unary_Transform(const Array1& a, TFuncT func, Array2& b) {
			thrust::transform(a.begin(), a.end(), b.begin(), func);
		}
		template<class TTFuncT, class Array1, class Array2, class Array3>
		void Binary_Transform(const Array1& a, const Array2& b, TTFuncT func, Array3& c) {
			thrust::transform(a.begin(), a.end(), b.begin(), c.begin(), func);
		}
	}

	namespace GPUFunc {
		template<class T> cudaDataType_t Cuda_Real_Type(void)
		{
			int siz = sizeof(T);
			if (siz == 4) { return CUDA_R_32F; }
			else if (siz == 8) { return CUDA_R_64F; }
			else { std::cerr << "[Error] AuxFuncCuda::Cuda_Type: Unknown data type\n"; return cudaDataType_t(); }
		}

	}

}