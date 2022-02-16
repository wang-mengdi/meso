#pragma once
#include "LinearMapping.h"
#include "Eigen/Dense"

namespace Meso {

	template<class T>
	class ConjugateGradient
	{
	public:
		LinearMapping<T>* linear_mapping = nullptr;//doesn't hold this value
		LinearMapping<T>* preconditioner = nullptr;//doesn't hold this value

		int max_iter = 0;
		real relative_tolerance = 1e-5;
		bool verbose = true;

		ArrayDv<T> b, x;

		//inner variables
		ArrayDv<T> p, Ap, z;


		cublasHandle_t cublasHandle = nullptr;

		ConjugateGradient() {	}
		~ConjugateGradient() {
			if (cublasHandle) cublasDestroy(cublasHandle);
		}

		//NOTE: it will take dof in Init() function and malloc accordingly
		void Init(LinearMapping<T>* _linear_mapping, LinearMapping<T>* preconditioner = nullptr, const int _max_iter = -1, const real _relative_tolerance = std::numeric_limits<real>::epsilon(), bool _verbose = false);

		////will reuse b to hold residule
		void Solve(ArrayDv<T>& x, ArrayDv<T>& b, int& iters, real& relative_error);
		//void Conjugate_Gradient(Scalar* x_dev, Scalar* b_dev, int& iters, Scalar& relative_error);
		////will reuse b_dev to store residual
		//void Solve_Device(Scalar* x_dev, Scalar* b_dev);
		//void Solve(Scalar* x_host, Scalar const* b_host);
	};

}