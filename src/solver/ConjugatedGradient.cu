//#include "ConjugatedGradient.h"
//#include <memory>
//#include <iostream>
//#include "cuda_runtime_api.h"
//#include "gpuUtils.h"
//
//void ConjugatedGradient::Init(const int _max_iter, const Scalar _relative_tolerance)
//{
//	Assert(linear_mapping != nullptr, "ConjugatedGradient::linear_mapping not initialized");
//	max_iter = _max_iter;
//	relative_tolerance = _relative_tolerance;
//	assert(linear_mapping->xDoF() == linear_mapping->yDoF());
//	if (preconditioner)
//	{
//		assert(preconditioner->xDoF() == preconditioner->yDoF());
//		assert(linear_mapping->xDoF() == preconditioner->xDoF());
//	}
//	dof = linear_mapping->xDoF();
//	cudaMalloc((void**)&b_dev, sizeof(Scalar) * dof);
//	cudaMalloc((void**)&x_dev, sizeof(Scalar) * dof);
//
//	cudaMalloc((void**)&d_p, sizeof(Scalar)*dof);
//	cudaMalloc((void**)&d_Ap, sizeof(Scalar)*dof);
//	cudaMalloc((void**)&d_z, sizeof(Scalar)*dof);
//
//	cublasCreate(&cublasHandle);
//}
//
//void ConjugatedGradient::Conjugate_Gradient(Scalar* d_x, Scalar* b_dev, int& iters, Scalar& relative_error)
//{
//	//https://flat2010.github.io/2018/10/26/%E5%85%B1%E8%BD%AD%E6%A2%AF%E5%BA%A6%E6%B3%95%E9%80%9A%E4%BF%97%E8%AE%B2%E4%B9%89/
//	//https://zhuanlan.zhihu.com/p/98642663
//	//See: simplex/docs/conjugate-gradient-notes-zh.md
//	Scalar one = 1, gamma;
//
//	//Use 0 as initial guess
//	//x0=0
//	cudaMemset(d_x, 0, sizeof(Scalar) * dof);
//	//initial residual is b
//	Scalar* d_r = b_dev;
//	
//	//int show_num = 100;
//	//Scalar* ptr_host = AuxFuncCPX::Global_Malloc<Scalar>(show_num, DataHolder::HOST);
//	//AuxFuncCPX::Global_Copy_Array(ptr_host, d_r, show_num, DataHolder::HOST, DataHolder::DEVICE);
//	//cudaDeviceSynchronize();
//	//std::cout << "value ptr extracted before solve: ";for (int i = 0;i < show_num;i++) std::cout << ptr_host[i] << " ";std::cout << "\n";
//	//AuxFuncCPX::Global_Free(ptr_host, DataHolder::HOST);
//
//
//	//rhs_norm2=r*r
//	Scalar rhs_norm2;
//	Dot(cublasHandle, dof, d_r, 1, d_r, 1, &rhs_norm2);
//	cudaDeviceSynchronize();
//	
//	//if b is zero, just solve to zero
//	if (rhs_norm2 == 0) {
//		//d_x is zero
//		iters = 0;
//		relative_error = 0;
//		return;
//	}
//	//(epsilon*|b|)^2
//	Scalar threshold_norm2 = relative_tolerance * relative_tolerance * rhs_norm2;
//
//	////z0=Minv*r0
//	if (preconditioner){preconditioner->applyMapping(d_z, d_r);}
//	else cudaMemcpy(d_z, d_r, sizeof(Scalar) * dof, cudaMemcpyDeviceToDevice);
//
//	//p0=z0
//	cudaMemcpy(d_p, d_z, sizeof(Scalar) * dof, cudaMemcpyDeviceToDevice);
//
//	//gamma0=dot(r0,z0)
//	Dot(cublasHandle, dof, d_z, 1, d_r, 1, &gamma);
//	cudaDeviceSynchronize();
//
//	Scalar fp;//fp_k=p_k^T*A*p_k
//	Scalar residual_norm2;//|r_k|^2
//	int i = 0;
//	for (i = 0;i < max_iter;i++) {
//		//Ap_k=A*p_k
//		linear_mapping->applyMapping(d_Ap, d_p);
//
//		//alpha_k=gamma_k/(p_k^T*A*p_k)
//		Dot(cublasHandle, dof, d_p, 1, d_Ap, 1, &fp);
//		cudaDeviceSynchronize();
//		Scalar alpha = gamma / fp;
//		Scalar neg_alpha = -alpha;
//
//		//x_{k+1} = x_k + alpha_k * p_k
//		//Axpy means y=y+a*x
//		Axpy(cublasHandle, dof, &alpha, d_p, 1, d_x, 1);
//
//		//r_{k+1} = r_k - alpha_k * Ap_k
//		Axpy(cublasHandle, dof, &neg_alpha, d_Ap, 1, d_r, 1);
//
//		Dot(cublasHandle, dof, d_r, 1, d_r, 1, &residual_norm2);
//		cudaDeviceSynchronize();
//		if (residual_norm2 < threshold_norm2) break;
//
//		//z_{k+1} = Minv * r_{k+1}
//		if (preconditioner) preconditioner->applyMapping(d_z, d_r);
//		else cudaMemcpy(d_z, d_r, sizeof(Scalar) * dof, cudaMemcpyDeviceToDevice);
//
//		//gamma_{k+1} = dot(r_{k+1}, z_{k+1})
//		Scalar gamma_old = gamma;
//		Dot(cublasHandle, dof, d_z, 1, d_r, 1, &gamma);
//		cudaDeviceSynchronize();
//		//if (gamma < threshold_norm2) break;
//
//		//beta_{k+1} = gamma_{k+1} / gamma_k
//		Scalar beta = gamma / gamma_old;
//
//		//p_{k+1} = z_{k+1} + beta_{k+1} * p_{k}
//		Scal(cublasHandle, dof, &beta, d_p, 1);
//		Axpy(cublasHandle, dof, &one, d_z, 1, d_p, 1);
//	}
//	iters = i;
//	relative_error = sqrt(residual_norm2 / rhs_norm2);
//}
//
//void ConjugatedGradient::Solve_Device(Scalar* x_dev, Scalar* b_dev)
//{
//	int iters;
//	Scalar relative_error;
//	Conjugate_Gradient(x_dev, b_dev, iters, relative_error);
//	if (verbose) Info("CPX CG take {} iterations with relative error {:.5e}", iters, relative_error);
//}
//
//void ConjugatedGradient::Solve(Scalar* x_host, Scalar const* b_host)
//{
//	int iters;
//	Scalar relative_error;
//	//AuxFuncCPX::Global_Copy_Array(b_dev, b_host, dof, DataHolder::DEVICE, DataHolder::HOST);
//	cudaMemcpy(b_dev, b_host, sizeof(Scalar) * dof, cudaMemcpyHostToDevice);
//	cudaDeviceSynchronize();
//
//	Solve_Device(x_dev, b_dev);
//
//	AuxFuncCPX::Global_Copy_Array(x_host, x_dev, dof, DataHolder::HOST, DataHolder::DEVICE);
//	cudaDeviceSynchronize();
//}
