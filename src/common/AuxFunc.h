//////////////////////////////////////////////////////////////////////////
// Common auxillary functions
// Copyright (c) (2022-), Bo Zhu, Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "Common.h"
#include "Constants.h"
#include <thrust/inner_product.h>
#include <thrust/transform.h>
#include <thrust/functional.h>
using namespace thrust::placeholders;

namespace Meso {
	namespace StringFunc{
		std::string To_String_Simple(const bool& a);

		template<class T,int d> std::string To_String_Simple(const Vector<T,d> & a) {
			std::string out = "[ ";
			if constexpr (std::is_same<T, C>::value == true) {
				for (int i = 0; i < d; i++) {
					C tmp_a(a[i]);
					out += "[";
					out += std::to_string(tmp_a.real());
					out += ", ";
					out += std::to_string(tmp_a.imag());
					out += "i]";
					if (i != d - 1) { out += " "; }
				}
				return out;
			}
			else {
				for (int i = 0; i < d; i++) {
					out += std::to_string(a[i]);
					if (i != d - 1) { out += " "; }
				}
			}
			out += ']';
			return out;
		}

		template<class T> std::string To_String_Simple(const T& a) {
			if constexpr (std::is_same<T, C>::value == true) {
				std::string out = "[";
				C tmp_a(a);
				out += std::to_string(tmp_a.real());
				out += ", ";
				out += std::to_string(tmp_a.imag());
				out += "i]";
				return out;
			}
			else {
				return std::to_string(a);
			}
		}

		Array<std::string> Split_String(const std::string& s, const std::string& delimiters = " \t\n\v\f\r");
	}

	namespace FileFunc {
		void Create_Directory(const bf::path path);
	}

	namespace OMPFunc {
		template<class IFunc>
		void Exec_Indices(const int n, IFunc f) {
#pragma omp parallel for
			for (int i = 0; i < n; i++) {
				f(i);
			}
		}
	}

	namespace MathFunc {
		template <typename T>  int Sign(T val) { return (T(0) < val) - (val < T(0)); }

		template<class T> __host__ __device__ T Eps(void) { return Eigen::NumTraits<T>::dummy_precision(); }
		template<class T> __host__ __device__ T Clamp(const T x, const T a, const T b) {
			if (x < a) return a;
			if (x > b) return b;
			return x;
		}
		////Scalar Math
		real Quick_Pow(real a, int n);//must ensure n>=0
		real Power2(const real a);//a^2
		real Power3(const real a);//a^3
		real Power4(const real a);//a^4
		real Power5(const real a);//a^5
		real Deg2Rad(real deg);

		////create vectors with compatible dimensions
		template<int d> __host__ __device__ Vector<real, d> V(const real x = (real)0, const real y = (real)0, const real z = (real)0);
		template<int d> __host__ __device__ Vector<int, d> Vi(const int x = 0, const int y = 0, const int z = 0, const int w = 0) {
			if constexpr (d == 1) { return Vector1i(x); }
			else if constexpr (d == 2) { return Vector2i(x, y); }
			else if constexpr (d == 3) { return Vector3i(x, y, z); }
			else if constexpr (d == 4) { return Vector4i(x, y, z, w); }
		}
		template<int d> __host__ __device__ Vector<real, d> V(const Vector2 v2);
		template<int d> __host__ __device__ Vector<real, d> V(const Vector3 v2);

		////Round up vector to a multiple of bn
		template<int d> Vector<int, d> Round_Up_To_Align(Vector<int, d> v, int bn) {
			for (int i = 0; i < d; i++) {
				//round to the nearest multiple of block_size
				v[i] = ((v[i] + bn - 1) / bn) * bn;
			}
			return v;
		}

		////Eigen zero compiler fix
		template<class T> __host__ __device__ constexpr T Zero() { return (T)0; }
#define Define_Zero(VecT) template<>  __host__ __device__ inline VecT Zero<VecT>(){return VecT::Zero();}
		Define_Zero(Matrix2); Define_Zero(Matrix3); Define_Zero(Matrix4);
		Define_Zero(Vector1); Define_Zero(Vector2); Define_Zero(Vector3);
		Define_Zero(Vector4);

		//atomic add two vectors
		template<class T, int d> __device__ __host__ void Atomic_Add(Vector<T,d>* a, Vector<T,d> b) {
			for (int i = 0; i < d; i++) { atomicAdd(&((*a)[i]), b[i]); }
		}

		//atomic add two values
		template<class T> __device__ __host__ void Atomic_Add(T* a, T b) {
			atomicAdd(a, b);
		}

		////Vector index sorting in ascending order
		template<int d> Vector<int, d> Sorted(const Vector<int, d>& _v) {
			Vector<int, d> v = _v;
			std::sort(v.data(), v.data() + v.size());
			return v;
		}

		Vector1 Orthogonal_Vector(const Vector1& v);
		Vector2 Orthogonal_Vector(const Vector2& v);	////this is the orthogonal vector on the *left* hand
		Vector3 Orthogonal_Vector(const Vector3& v);
		Vector4 Orthogonal_Vector(const Vector4& v);

		Vector1 Cross(const Vector2& v1, const Vector2& v2);
		Vector2 Cross(const Vector1& v1, const Vector2& v2);

		real Angle_Between(const Vector2& v1, const Vector2& v2);	////[0,pi]
		real Angle_Between(const Vector3& v1, const Vector3& v2);	////[0,pi]
		real Angle_From_To(const Vector2& v1, const Vector2& v2);	////[-pi,pi], default is +z
		real Angle_From_To(const Vector3& v1, const Vector3& v2);	////[-pi,pi], no axis specified, [0,pi]

		//p0-p1 is the shared edge, range: [-pi,+pi]
		inline real Dihedral_Angle(const Vector3& n0, const Vector3& n1, const Vector3& p0, const Vector3& p1) {
			real cosine = n0.dot(n1);
			Vector3 e = (p0 - p1).normalized();
			real sine = e.dot(n0.cross(n1));
			return atan2(sine, cosine);
		}

		inline Vector2 Barycentric_Weights(const Vector3& p, const Vector3& p0, const Vector3& p1) {
			Vector3 e = p1 - p0;
			double t = e.dot(p - p0) / e.dot(e);
			return Vector2(1 - t, t);
		}

		//distance from one vertex to base, p0-p1 is the base, always positive
		inline real Distance(const Vector3& p, const  Vector3& p0, const  Vector3& p1) {
			return ((p - p0).cross(p - p1)).norm() / (p0 - p1).norm();
		}

		template<class T, int dim> inline T Abs_Min(const Vector<T, dim>& v) { T v_min = abs(v[0]); for (int i = 1; i < dim; i++)if (abs(v[i]) < v_min)v_min = v[i]; return v_min; }
		template<class T, int dim> inline T Abs_Max(const Vector<T, dim>& v) { T v_max = abs(v[0]); for (int i = 1; i < dim; i++)if (abs(v[i]) > v_max)v_max = v[i]; return v_max; }
		template<class T, int dim> int Abs_Min_Index(const Vector<T, dim>& v) { int i_min = 0; for (int i = 1; i < dim; i++)if (abs(v[i]) < abs(v[i_min]))i_min = i; return i_min; }
		template<class T, int dim> int Abs_Max_Index(const Vector<T, dim>& v) { int i_max = 0; for (int i = 1; i < dim; i++)if (abs(v[i]) > abs(v[i_max]))i_max = i; return i_max; }
		template<class T, int dim> inline int Min_Index(const Vector<T, dim>& v) { int i_min = 0; for (int i = 1; i < dim; i++)if (v[i] < v[i_min])i_min = i; return i_min; }
		template<class T, int dim> inline int Max_Index(const Vector<T, dim>& v) { int i_max = 0; for (int i = 1; i < dim; i++)if (v[i] > v[i_max])i_max = i; return i_max; }
		template<class T, int dim> inline Vector<T, dim> Cwise_Min(const Vector<T, dim>& v1, const Vector<T, dim>& v2) { return v1.cwiseMin(v2); }
		template<class T, int dim> inline Vector<T, dim> Cwise_Max(const Vector<T, dim>& v1, const Vector<T, dim>& v2) { return v1.cwiseMax(v2); }
		template<class T, int dim> inline Vector<T, dim> Cwise_Multiply(const Vector<T, dim>& v1, const Vector<T, dim>& v2) { return v1.cwiseProduct(v2); }
		template<class T, int dim> inline Vector<T, dim> Cwise_Divide(const Vector<T, dim>& v1, const Vector<T, dim>& v2) { Vector<T, dim> v; for (int i = 0; i < dim; i++) v[i] = v1[i] / v2[i]; return v; }

		////eigenvectors, svd, polar decomposition
		////eigenvector/value corresponding to the eigenvalue with the max magnitude
		Vector2 Principal_Eigenvector(const Matrix2& v);
		Vector3 Principal_Eigenvector(const Matrix3& v);
		real Principal_Eigenvalue(const Matrix2& v);
		real Principal_Eigenvalue(const Matrix3& v);

		////eigenvector/value corresponding to the eigenvalue with the min magnitude
		Vector2 Min_Eigenvector(const Matrix2& v);
		Vector3 Min_Eigenvector(const Matrix3& v);
		real Min_Eigenvalue(const Matrix2& v);
		real Min_Eigenvalue(const Matrix3& v);

		//absolute(a - b) <= (atol + rtol * absolute(b))
		template<typename DerivedA, typename DerivedB>
		bool All_Close(const Eigen::DenseBase<DerivedA>& a,
			const Eigen::DenseBase<DerivedB>& b,
			const typename DerivedA::RealScalar& rtol
			= Eigen::NumTraits<typename DerivedA::RealScalar>::dummy_precision(),
			const typename DerivedA::RealScalar& atol
			= Eigen::NumTraits<typename DerivedA::RealScalar>::epsilon())
		{
			return ((a.derived() - b.derived()).array().abs()
				<= (atol + rtol * b.derived().array().abs())).all();
		}

		template<typename T1, typename T2> bool Close(T1 a, T2 b, const T1 rtol
			= Eigen::NumTraits<T1>::dummy_precision(),
			const T1 atol
			= Eigen::NumTraits<T1>::epsilon()) {
			return abs(a - b) <= atol + rtol * abs(b);
		}
	}

	namespace ArrayFunc {
		//Comparison，array or vector
		template<class ArrayT> bool All_Less(const ArrayT& a0, const ArrayT& a1) { for (auto i = 0; i < a0.size(); i++) { if (a0[i] >= a1[i])return false; }return true; }
		template<class ArrayT> bool All_Less_Equal(const ArrayT& a0, const ArrayT& a1) { for (auto i = 0; i < a0.size(); i++) { if (a0[i] > a1[i])return false; }return true; }
		template<class ArrayT> bool All_Greater(const ArrayT& a0, const ArrayT& a1) { for (auto i = 0; i < a0.size(); i++) { if (a0[i] <= a1[i])return false; }return true; }
		template<class ArrayT> bool All_Greater_Equal(const ArrayT& a0, const ArrayT& a1) { for (auto i = 0; i < a0.size(); i++) { if (a0[i] < a1[i])return false; }return true; }
		template<class ArrayT> bool Has_Equal(const ArrayT& a0, const ArrayT& a1) { for (auto i = 0; i < a0.size(); i++)if (a0[i] == a1[i])return true; return false; }
		template<class ArrayT> bool Has_Less_Equal(const ArrayT& a0, const ArrayT& a1) { for (auto i = 0; i < a0.size(); i++)if (a0[i] <= a1[i])return true; return false; }
		template<class ArrayT> bool Has_Greater_Equal(const ArrayT& a0, const ArrayT& a1) { for (auto i = 0; i < a0.size(); i++)if (a0[i] >= a1[i])return true; return false; }

		template<class Array1>
		constexpr auto Data(Array1& arr)noexcept {
			return thrust::raw_pointer_cast(arr.data());
		}

		template<class T, DataHolder side>
		bool Has_Zero(const Array<T, side>& a) {
			if constexpr (side == HOST) {
				for (int i = 0; i < a.size(); i++) if (a[i] == 0) return true;
			}
			else {
				Array<T, HOST> b = a;
				for (int i = 0; i < b.size(); i++) if (b[i] == 0) return true;
			}
			return false;
		}
		template<class T, DataHolder side>
		bool Is_Finite(const Array<T, side>& a) {
			if constexpr (side == HOST) {
				for (int i = 0; i < a.size(); i++) if (!std::isfinite(a[i])) return false;
			}
			else {
				Array<T, HOST> b = a;
				for (int i = 0; i < b.size(); i++) if (!std::isfinite(b[i])) return false;
			}
			return true;
		}
		template<class T>
		bool Equals(const Array<T, HOST>& a, decltype(a) b) {
			if (a.size() != b.size()) return false;
			for (int i = 0; i < a.size(); i++) if (a[i] != b[i]) return false;
			return true;
		}
		template<class T, DataHolder side = HOST>
		T Max_Abs(const Array<T, side>& a) {
			return thrust::reduce(
				a.begin(),
				a.end(),
				(T)0,
				[=](const T a, const T b) { return std::max(std::abs(a), std::abs(b)); }
			);
		}
		template<class Array1, class T>
		void Fill(Array1& a, const T val) {
			thrust::fill(a.begin(), a.end(), val);
		}
		template<class Array1>
		double Dot(const Array1& a, decltype(a) b) {
			Assert(a.size() == b.size(), "[ArrayFunc::Dot] try to dot length {} against {}", a.size(), b.size());
			return thrust::inner_product(a.begin(), a.end(), b.begin(), (double)0);
		}

		template<class Array1> __host__ __device__
			C Conj_Dot(const Array1& a, const Array1& b) {
			//printf("hit conj dot!\n");
			//Assert(a.size() == b.size(), "[ArrayFunc::Conj_Dot] try to dot length {} against {}", a.size(), b.size());
			//printf("hit after assert!\n");
			//printf("a size: %d b size: %d\n", a.size(), b.size());
			struct a_dot_conj_b : public thrust::binary_function<C, C, C>
			{
				__host__ __device__
					C operator()(C a, C b)
				{
					return  a * thrust::conj(b);
				};
			};
			//printf("hit conj dot!\n");
			C ans = thrust::inner_product(a.begin(), a.end(), b.begin(), C(0.0, 0.0), thrust::plus<C>(), a_dot_conj_b());
			//printf("ans: (%f,%f)\n",ans.real(), ans.imag());
			return ans;

		}
		template<class Array1>
		double Norm(const Array1& a) {
			double squared_norm = Dot(a, a);
			return sqrt(squared_norm);
		}
		template<class T>
		int Largest_Norm_Element(const Array<T>& arr)
		{
			int idx = -1; real max_norm = -1.0;
			for (int i = 0; i < arr.size(); i++) {
				real v = arr[i].norm();
				if (v >= max_norm) {
					idx = i;
					max_norm = v;
				}
			}
			return idx;
		}
		template<class T> 
		real Largest_Norm(const Array<T>& arr) {
			int idx = Largest_Norm_Element<T>(arr);
			if (idx < 0) return 0;
			else return arr[idx].norm();
		}
		template<class Array1>
		decltype(auto) Sum(const Array1& v) {
			using T = decltype(v[0]);
			return thrust::reduce(
					v.begin(),
					v.end(),
					(T)0,
					thrust::plus<T>()
				);
		}
		template<class T, DataHolder side = HOST>
		T Mean(const Array<T, side>& v) { 
			return 1. / (real)v.size() *
				thrust::reduce(
					v.begin(),
					v.end(),
					(T)0,
					thrust::plus<T>()
				);
		}
		template<class T> 
		T Abs_Mean(const Array<T>& v) { 
			return 1. / (real)v.size() *
				thrust::reduce(
					v.begin(),
					v.end(),
					(T)0,
					[=](const T a, const T b) { return std::abs(a) + std::abs(b); }
			);
		}
		//scalar multiplication, a*=b (where b is a scalar)
		template<class Array1, class T>
		void Multiply_Scalar(Array1& a, const T b) {
			thrust::transform(a.begin(), a.end(), a.begin(), _1 * b);
		}
		//element-wise multiplication, a*.=b
		template<class Array1, class Array2>
		void Multiply(Array1& a, const Array2 &b) {
			thrust::transform(a.begin(), a.end(), b.begin(), a.begin(), _1 * _2);
		}
		//element-wise divide, a/.=b
		template<class Array1, class Array2>
		void Divide(Array1& a, const Array2& b) {
			thrust::transform(a.begin(), a.end(), b.begin(), a.begin(), _1 / _2);
		}
		//add a scalar
		template<class Array1, class T>
		void Add_Scalar(Array1& a, const T b) { thrust::transform(a.begin(), a.end(), a.begin(), _1 + b); }
		//element-wise add, a+.=b
		template<class Array1, class Array2>
		void Add(Array1& a, const Array2& b) { thrust::transform(a.begin(), a.end(), b.begin(), a.begin(), _1 + _2); }
		template<class Array1, class T>
		void operator += (Array1& a, const T b) {
			thrust::transform(a.begin(), a.end(), a.begin(), _1 + b);
		}
		template<class Array1, class Array2>
		void Minus(Array1& a, const Array2& b) {
			thrust::transform(a.begin(), a.end(), b.begin(), a.begin(), _1 - _2);
		}
		template<class Array1, class Array2, class T>
		void Set_Masked(Array1& a, const Array2& b, const T val) {
			thrust::transform_if(
				a.begin(),//first
				a.end(),//last
				b.begin(),//stencil
				a.begin(),//result
				[=]__host__ __device__(const T v1)->T {return val; },//op
				thrust::identity<bool>()//pred
			);
		}
		//a=b, note it's reverse order of thrust::copy itself
		template<class Array1, class Array2>
		void Copy(Array1& a, const Array2& b) {
			thrust::copy(b.begin(), b.end(), a.begin());
		}
		//copy all places with mask==true
		template<class Array1, class Array2, class Array3>
		void Copy_Masked(Array1& dest, const Array2& src, const Array3& mask) {
			thrust::transform_if(
				src.begin(),//first
				src.end(),//last
				mask.begin(),//stencil
				dest.begin(),//result
				_1,//op
				thrust::identity<bool>()//pred
			);
		}
		//copy all places with mask==false
		template<class Array1, class Array2, class Array3>
		void Copy_UnMasked(Array1& dest, const Array2& src, const Array3& mask) {
			thrust::transform_if(
				src.begin(),//first
				src.end(),//last
				mask.begin(),//stencil
				dest.begin(),//result
				_1,//op
				thrust::logical_not<bool>()//pred
			);
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
		template<class T>
		bool Is_Approx(const Array<T>& a, const Array<T>& b) {
			//similar to isApprox() in Eigen
			if (a.size() != b.size()) return false;
			double a_norm2 = Dot(a, a);
			double b_norm2 = Dot(b, b);
			Array<T> res; res = a;
			Minus<Array<T>>(res, b);
			double res_norm2 = Dot(res, res);
			double p = Eigen::NumTraits<T>::dummy_precision();
			return res_norm2 <= p * p * std::min(a_norm2, b_norm2);
		}
		template<class T>
		bool Is_Approx(const ArrayDv<T>& a, const ArrayDv<T>& b) {
			Array<T> a_host = a, b_host = b;
			return Is_Approx<T>(a_host, b_host);
		}
//		////Dim conversion for vectors and vector arrays
//		template<class T, int d1, int d2> void Dim_Conversion(const Vector<T, d1>& input, Vector<T, d2>& output, const T filled_value = (T)0)
//		{
//			constexpr int n = d1 < d2 ? d1 : d2;
//			for (int i = 0; i < n; i++)output[i] = input[i];
//			if /*constexpr*/ (n < d2) {
//				for (int i = n; i < d2; i++)output[i] = filled_value;
//			}
//		}
//
//		template<class T, int d1, int d2> void Dim_Conversion_Array(const Array<Vector<T, d1> >& input, Array<Vector<T, d2> >& output, const T filled_value = (T)0)
//		{
//			const int n = (int)input.size();
//#pragma omp parallel for
//			for (auto i = 0; i < n; i++) { Dim_Conversion<T, d1, d2>(input[i], output[i], filled_value); }
//		}
//
//		////Dim conversion for matrices and matrix arrays
//		template<class T, int d1, int d2> void Dim_Conversion(const Matrix<T, d1>& input, Matrix<T, d2>& output, const T filled_value = (T)0)
//		{
//			constexpr int n = d1 < d2 ? d1 : d2;
//			output = Matrix<T, d2>::Constant(filled_value);
//			for (int i = 0; i < n; i++)for (int j = 0; j < n; j++)output(i, j) = input(i, j);
//		}
//
//		template<class T, int d1, int d2> void Dim_Conversion_Array(const Array<Matrix<T, d1> >& input, Array<Matrix<T, d2> >& output)
//		{
//			const int n = (int)input.size();
//#pragma omp parallel for
//			for (auto i = 0; i < n; i++)Dim_Conversion<T, d1, d2>(input[i], output[i]);
//		}
	}

	namespace GPUFunc {
		template<class T> cudaDataType_t Cuda_Real_Type(void)
		{
			int siz = sizeof(T);
			if (siz == 4) { return CUDA_R_32F; }
			else if (siz == 8) { return CUDA_R_64F; }
			else { std::cerr << "[Error] AuxFuncCuda::Cuda_Type: Unknown data type\n"; return cudaDataType_t(); }
		}

		template<int d>
		__device__ __host__ Vector<int, d> Thread_Coord(const dim3 blockIdx, const dim3 threadIdx) {
			if constexpr (d == 2) {
				return Vector2i(blockIdx.x * 8 + threadIdx.x, blockIdx.y * 8 + threadIdx.y);
			}
			else if constexpr (d == 3) {
				return Vector3i(blockIdx.x * 4 + threadIdx.x, blockIdx.y * 4 + threadIdx.y, blockIdx.z * 4 + threadIdx.z);
			}
		}

		template<typename A, typename F>
		__global__ void Cwise_Mapping(A v1, F f, int N)
		{
			int i = blockIdx.x * blockDim.x + threadIdx.x;
			if (i >= N) return;
			f(v1[i]);
		}

		template<typename A, typename F>
		void Cwise_Mapping_Wrapper(A v1, F f, int N)
		{
			Cwise_Mapping << <((N + 63) >> 6), 64 >> > (v1, f, N);
		}

		template<typename A, typename B, typename F>
		__global__ void Cwise_Mapping(A v1, B v2, F f, int N)
		{
			int i = blockIdx.x * blockDim.x + threadIdx.x;
			if (i >= N) return;
			f(v1[i], v2[i]);
		}

		template<typename A, typename B, typename F>
		void Cwise_Mapping_Wrapper(A v1, B v2, F f, int N)
		{
			Cwise_Mapping << <((N + 63) >> 6), 64 >> > (v1, v2, f, N);
		}

		template<typename A, typename B, typename C, typename F>
		__global__ void Cwise_Mapping(A v1, B v2, C v3, F f, int N)
		{
			int i = blockIdx.x * blockDim.x + threadIdx.x;
			if (i >= N) return;
			f(v1[i], v2[i], v3[i]);
		}

		template<typename A, typename B, typename C, typename F>
		void Cwise_Mapping_Wrapper(A v1, B v2, C v3, F f, int N)
		{
			Cwise_Mapping << <((N + 63) >> 6), 64 >> > (v1, v2, v3, f, N);
		}
	}

}