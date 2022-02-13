//#pragma once
//#include "LinearMapping.h"
//#include "Eigen/Dense"
//#include "cublas_v2.h"
//#include "AuxFuncCPX.h"
//
//using namespace AuxFuncCPX;
//
//class ConjugatedGradient
//{
//public:
//	LinearMapping* linear_mapping = nullptr;//doesn't hold this value
//	LinearMapping *preconditioner = nullptr;//doesn't hold this value
//	int dof = 0;
//	Scalar* b_dev = nullptr, *x_dev = nullptr;
//	//inner variables
//	Scalar* d_p = nullptr, * d_Ap = nullptr, * d_z = nullptr;
//
//	Scalar relative_tolerance = (Scalar)1e-5;
//
//	bool verbose = true;
//	int max_iter = 0;
//
//	cublasHandle_t cublasHandle = nullptr;
//
//	ConjugatedGradient() {	}
//	~ConjugatedGradient() {
//		//note that Global_Free() can automatically handle nullptr
//		Global_Free(b_dev, DataHolder::DEVICE);
//		Global_Free(x_dev, DataHolder::DEVICE);
//
//		Global_Free(d_p, DataHolder::DEVICE);
//		Global_Free(d_Ap, DataHolder::DEVICE);
//		Global_Free(d_z, DataHolder::DEVICE);
//		if (cublasHandle) cublasDestroy(cublasHandle);
//	}
//
//	//NOTE: it will take dof in Init() function and malloc accordingly
//	void Init(const int _max_iter, const Scalar _relative_tolerance = 1e-5);
//
//	//will reuse b_dev to store residual
//	void Conjugate_Gradient(Scalar* x_dev, Scalar* b_dev, int& iters, Scalar& relative_error);
//	//will reuse b_dev to store residual
//	void Solve_Device(Scalar* x_dev, Scalar* b_dev);
//	void Solve(Scalar* x_host, Scalar const* b_host);
//};