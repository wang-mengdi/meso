#pragma once
#include "LinearMapping.h"
#include "Eigen/Dense"
#include "cublas_v2.h"

template<class T>
class ConjugatedGradient
{
public:
	LinearMapping<T>* linear_mapping = nullptr;//doesn't hold this value
	LinearMapping<T>* preconditioner = nullptr;//doesn't hold this value
	
	int max_iter = 0;
	T relative_tolerance = 1e-5;
	bool verbose = true;

	ArrayDv<T> b, x;
	
	//inner variables
	int dof;
	ArrayDv<T> p, Ap, z;
	

	cublasHandle_t cublasHandle = nullptr;

	ConjugatedGradient() {	}
	~ConjugatedGradient() {
		if (cublasHandle) cublasDestroy(cublasHandle);
	}

	//NOTE: it will take dof in Init() function and malloc accordingly
	void Init(LinearMapping<T>* _linear_mapping, LinearMapping<T>* preconditioner = nullptr, const int _max_iter = -1, const T _relative_tolerance = std::numeric_limits<T>::epsilon(), bool _verbose = false);

	////will reuse b_dev to store residual
	//void Conjugate_Gradient(Scalar* x_dev, Scalar* b_dev, int& iters, Scalar& relative_error);
	////will reuse b_dev to store residual
	//void Solve_Device(Scalar* x_dev, Scalar* b_dev);
	//void Solve(Scalar* x_host, Scalar const* b_host);
};