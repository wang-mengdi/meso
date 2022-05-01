//////////////////////////////////////////////////////////////////////////
// Sparse solver
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __Sparse_h__
#define __Sparse_h__
#include <Eigen/Sparse>
#include "Common.h"

using namespace Meso;
////Eigen sparse type alias
template<class T> using VectorN=Eigen::Matrix<T,-1,1>;
template<class T> using InnerIterator=typename SparseMatrix<T>::InnerIterator;
template<class T> using DiagonalMatrix=Eigen::DiagonalMatrix<T,Eigen::Dynamic,Eigen::Dynamic>;
template<class T> using Triplet=Eigen::Triplet<T,int>;

namespace SparseFunc{
////block matrix operations
template<int dim,class T_MAT> void Add_Block(SparseMatrix<real>& K,const int K_i,const int K_j,const T_MAT& K_b,const int Kb_i=0,const int Kb_j=0)
{for(int i=0;i<dim;i++)for(int j=0;j<dim;j++){K.coeffRef(K_i*dim+i,K_j*dim+j)+=K_b.coeff(Kb_i*dim+i,Kb_j*dim+j);}}

template<int dim,class T_MAT> void Copy_Block(SparseMatrix<real>& K,const int K_i,const int K_j,const T_MAT& K_b,const int Kb_i=0,const int Kb_j=0)
{for(int i=0;i<dim;i++)for(int j=0;j<dim;j++){K.coeffRef(K_i*dim+i,K_j*dim+j)=K_b.coeff(Kb_i*dim+i,Kb_j*dim+j);}}

template<int dim> void Set_Block(SparseMatrix<real>& K,const int K_i,const int K_j,const real value)
{for(int i=0;i<dim;i++)for(int j=0;j<dim;j++){K.coeffRef(K_i*dim+i,K_j*dim+j)=value;}}

inline void Set_Value(SparseMatrix<real>& K,const real value)	/////set all nonzero values to be value but keep all indices
{
	int nnz=(int)K.nonZeros();
	#pragma omp parallel for
	for(int i=0;i<nnz;i++){K.valuePtr()[i]=value;}
}
};

////default Eigen sparse (and dense) linear solvers
namespace SparseSolver{
    class Params
    {public:
        real tolerance=(real)1e-5;
        int max_iter_num=1000;
    };

    bool Conjugate_Gradient(const SparseMatrix<real>& A,VectorN<real>& x,const VectorN<real>& b,const Params params=Params());
	bool IC_PCG(const SparseMatrix<real>& A,VectorN<real>& x,const VectorN<real>& b,const Params params=Params());
    bool BiCG_STAB(const SparseMatrix<real>& A,VectorN<real>& x,const VectorN<real>& b,const Params params=Params());
    template<class T> bool LU(const SparseMatrix<T>& A,VectorN<T>& x,const VectorN<T>& b,const Params params=Params());
	template<class T> bool Least_Squares(const SparseMatrix<T>& A,VectorN<T>& x,const VectorN<T>& b,const Params params=Params());
};
#endif