//////////////////////////////////////////////////////////////////////////
// Random number
// Copyright (c) (2022-), Bo Zhu, Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once

#include <random>
#include <chrono>
#include "Common.h"

namespace Meso {

	//We will imitate interfaces of Python's random lib here
	namespace Random {
		//int-valued functions
		int RandInt(int a, int b);//random number in [a,b]

		//real-valued functions
		real Random(void);//random number in [0,1)
		real Uniform(real a, real b);

        template<class T>
        Array<T> Random_Array(const int n, const T a = 0.0, const decltype(a) b = 1.0) {
            Array<T> arr(n);
            for (int i = 0; i < n; i++) arr[i] = Uniform(a, b);
            return arr;
        }
        VectorXd Random_VectorXd(int n, real a = 0.0, real b = 1.0);

        template<class T>
        void Sparse_Diagonal_Dominant_Matrix(int cols, int rows, Eigen::SparseMatrix<T, Eigen::RowMajor, int>& mat) {
            std::vector<Eigen::Triplet<T>> tripletList;
            tripletList.reserve(3 * rows); //three entries per row
            for (int i = 0; i < rows; i++)
            {
                for (int j = std::max(i - 1, 0); j <= std::min(i + 1, cols - 1); j++) {
                    if (i == j) { tripletList.push_back(Eigen::Triplet<T>(i, j, Random::Random() * 100)); }
                    else { tripletList.push_back(Eigen::Triplet<T>(i, j, Random::Random())); }
                }
            }
            mat.resize(rows, cols);
            mat.setFromTriplets(tripletList.begin(), tripletList.end());
        }
	}

}