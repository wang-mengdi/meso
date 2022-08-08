#pragma once
#include "LinearMapping.h"
#include "Eigen/Dense"
#include <tuple>

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
		bool is_pure_neumann = false;
		//inner variables
		ArrayDv<T> p, Ap, z;


		cublasHandle_t cublasHandle = nullptr;

		ConjugateGradient() {	}
		~ConjugateGradient() {
			if (cublasHandle) cublasDestroy(cublasHandle);
		}

		//NOTE: it will take dof in Init() function and malloc accordingly
		void Init(LinearMapping<T>* _linear_mapping, LinearMapping<T>* _preconditioner = nullptr, bool _verbose = false, const int _max_iter = -1, const real _relative_tolerance = std::numeric_limits<real>::epsilon()) {
			linear_mapping = _linear_mapping;
			preconditioner = _preconditioner;
			Assert(linear_mapping != nullptr, "[ConjugateGradient] linear_mapping not initialized");
			if (_max_iter == -1) max_iter = linear_mapping->XDoF() * 2;
			else max_iter = _max_iter;
			relative_tolerance = _relative_tolerance;
			verbose = _verbose;

			Assert(linear_mapping->XDoF() == linear_mapping->YDoF(), "[ConjugateGradient] row number and col number must be equal");
			if (preconditioner)	Assert(
				preconditioner->XDoF() == preconditioner->YDoF() && linear_mapping->XDoF() == preconditioner->XDoF(),
				"[ConjugateGradient] preconditioner size must be equal to matrix size"
			);

			int dof = linear_mapping->XDoF();
			p.resize(dof);
			Ap.resize(dof);
			z.resize(dof);

			if (cublasHandle) cublasDestroy(cublasHandle);
			cublasCreate(&cublasHandle);
		}

		////return (iters,relative_error)
		////will destroy b and reuse it to hold residule
		std::tuple<int, real> Solve(ArrayDv<T>& x, ArrayDv<T>& b) {
			//https://flat2010.github.io/2018/10/26/%E5%85%B1%E8%BD%AD%E6%A2%AF%E5%BA%A6%E6%B3%95%E9%80%9A%E4%BF%97%E8%AE%B2%E4%B9%89/
			//https://zhuanlan.zhihu.com/p/98642663
			//See: docs/mgpcg-notes-zh.md

			//Use 0 as initial guess
			//x0=0
			x.resize(linear_mapping->XDoF());
			thrust::fill(x.begin(), x.end(), 0);
			//initial residual is b

			auto& r = b;

			//treatment for pure neumann 
			if (is_pure_neumann)
			{
				T r_mean = ArrayFunc::Mean<T, DEVICE>(r);
				ArrayDv<T> mu;
				mu.resize(linear_mapping->XDoF());
				thrust::fill(mu.begin(), mu.end(), r_mean);
				ArrayFunc::Minus(r, mu);
			}

			//rhs_norm2=r*r
			double rhs_norm2 = ArrayFunc::Dot(r, r);
			if (verbose) Info("ConjugateGradient initial norm of rhs: {}", sqrt(rhs_norm2));

			//if b is zero, just solve to zero
			if (rhs_norm2 == 0) {
				//d_x is zero
				//iters=0, relative_error=0
				return std::make_tuple(0, (real)0);
			}
			//(epsilon*|b|)^2
			double threshold_norm2 = relative_tolerance * relative_tolerance * rhs_norm2;
			threshold_norm2 = std::max(threshold_norm2, std::numeric_limits<double>::min());


			////z0=Minv*r0
			if (preconditioner) preconditioner->Apply(z, r);
			else ArrayFunc::Copy(z, r);

			//p0=z0
			ArrayFunc::Copy(p, z);

			//gamma0=dot(r0,z0)
			double gamma = ArrayFunc::Dot(z, r);

			double residual_norm2;//|r_k|^2
			int i = 0;
			for (i = 0; i < max_iter; i++) {
				//Ap_k=A*p_k
				linear_mapping->Apply(Ap, p);

				//alpha_k=gamma_k/(p_k^T*A*p_k)
				double fp = ArrayFunc::Dot(p, Ap);//fp_k=p_k^T*A*p_k
				double alpha = gamma / fp;

				Assert(std::isnormal(alpha), "ConjugateGradient: alpha={} at iter {}", alpha, i);

				//Info("iter {} alpha {}", i, alpha);

				//x_{k+1} = x_k + alpha_k * p_k
				//Axpy means y=y+a*x
				ArrayFunc::Axpy(alpha, p, x);

				//r_{k+1} = r_k - alpha_k * Ap_k
				ArrayFunc::Axpy(-alpha, Ap, r);

				//treatment for pure neumann 
				if (is_pure_neumann)
				{
					T r_mean = ArrayFunc::Mean<T, DEVICE>(r);
					ArrayDv<T> mu;
					mu.resize(linear_mapping->XDoF());
					thrust::fill(mu.begin(), mu.end(), r_mean);
					ArrayFunc::Minus(r, mu);
				}

				residual_norm2 = ArrayFunc::Dot(r, r);
				if (verbose) Info("ConjugateGradient iter {} norm {} against threshold {}", i, sqrt(residual_norm2), sqrt(threshold_norm2));
				if (residual_norm2 < threshold_norm2) break;

				//z_{k+1} = Minv * r_{k+1}
				if (preconditioner) preconditioner->Apply(z, r);
				else ArrayFunc::Copy(z, r);

				//gamma_{k+1} = dot(r_{k+1}, z_{k+1})
				double gamma_old = gamma;
				gamma = ArrayFunc::Dot(z, r);

				//beta_{k+1} = gamma_{k+1} / gamma_k
				double beta = gamma / gamma_old;

				//p_{k+1} = z_{k+1} + beta_{k+1} * p_{k}
				ArrayFunc::Scal(beta, p);
				ArrayFunc::Axpy(1, z, p);
			}

			//return (iters,relative_error)
			return std::make_tuple(i, (real)sqrt(residual_norm2 / rhs_norm2));
		}
	};

}