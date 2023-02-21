#pragma once
#include "PoissonMapping.h"
#include "Eigen/Dense"
#include <tuple>

namespace Meso {

	template<class T, int d>
	class ConjugateGradient
	{
	public:
		MaskedPoissonMapping<T, d>* poisson_ptr = nullptr;//doesn't hold this value
		LinearMapping<T>* preconditioner = nullptr;//doesn't hold this value

		int max_iter = 0;
		real relative_tolerance = 1e-5;
		bool verbose = true;
		bool pure_neumann = false;
		int dof;
		//inner variables
		ArrayDv<T> p, Ap, z, mu;


		cublasHandle_t cublasHandle = nullptr;

		ConjugateGradient() {	}
		~ConjugateGradient() {
			if (cublasHandle) cublasDestroy(cublasHandle);
		}

		//NOTE: it will take dof in Init() function and malloc accordingly
		void Init(MaskedPoissonMapping<T, d>* _poisson_ptr, LinearMapping<T>* _preconditioner = nullptr, bool _verbose = false, const int _max_iter = -1, const real _relative_tolerance = std::numeric_limits<T>::epsilon(), bool _pure_neumann = false) {
			poisson_ptr = _poisson_ptr;
			preconditioner = _preconditioner;
			Assert(poisson_ptr != nullptr, "[ConjugateGradient] poisson_ptr not initialized");
			if (_max_iter == -1) max_iter = poisson_ptr->XDoF() * 2;
			else max_iter = _max_iter;
			relative_tolerance = _relative_tolerance;
			verbose = _verbose;
			pure_neumann = _pure_neumann;

			Assert(poisson_ptr->XDoF() == poisson_ptr->YDoF(),"[ConjugateGradient] row number and col number must be equal");
			if (preconditioner)	Assert(
				preconditioner->XDoF() == preconditioner->YDoF() && poisson_ptr->XDoF() == preconditioner->XDoF(),
				"[ConjugateGradient] preconditioner size must be equal to matrix size"
			);

			dof = poisson_ptr->XDoF();
			p.resize(dof);
			Ap.resize(dof);
			z.resize(dof);
			if (pure_neumann)
				mu.resize(dof);

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
			x.resize(dof);
			thrust::fill(x.begin(), x.end(), 0);
			//initial residual is b

			auto& r = b;

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

			if (pure_neumann)
				Minus_Average(r, poisson_ptr->cell_type, mu, dof);

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
				poisson_ptr->Apply(Ap, p);
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

				residual_norm2 = ArrayFunc::Dot(r, r);

				if (verbose) Info("ConjugateGradient iter {} norm {} against threshold {}", i, sqrt(residual_norm2), sqrt(threshold_norm2));
				if (residual_norm2 < threshold_norm2) break;

				if (pure_neumann)
					Minus_Average(r, poisson_ptr->cell_type, mu, dof);

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

		void Minus_Average(ArrayDv<T>& _r, const FieldDv<unsigned char, d>& _cell_type, ArrayDv<T>& _mu, const int _dof)
		{
			
			T sum = ArrayFunc::Sum<T, DEVICE>(_r);
			int fluid_cnt = ArrayFunc::Count<unsigned char, DEVICE>(*(_cell_type.data), 0);
			T val = sum / fluid_cnt;
			auto cond_set = [val]__device__(T & tv1, const unsigned char type) { if (type == 0) tv1 = val; else tv1 = 0; };
			T* mu_ptr = thrust::raw_pointer_cast(_mu.data());
			const unsigned char* cell_type_ptr = _cell_type.Data_Ptr();
			GPUFunc::Cwise_Mapping_Wrapper(mu_ptr, cell_type_ptr, cond_set, _dof);
			ArrayFunc::Minus(_r, _mu);
		}
	};
}