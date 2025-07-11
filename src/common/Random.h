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
        template<class T, int d>
        Vector<T, d> Uniform_In_Box(const Vector<T, d> a, decltype(a) b) {
            Vector<T, d> ret;
            for (int i = 0; i < d; i++) ret[i] = Uniform(a[i], b[i]);
            return ret;
        }

        template<class T>
        void Fill_Random_Array(Array<T>& arr, const T a = 0.0, const T b = 1.0) {
            for (int i = 0; i < arr.size(); i++) {
                arr[i] = Uniform(a, b);
            }
        }
        template<class T>
        Array<T> Random_Array(const int n, const T a = 0.0, const decltype(a) b = 1.0) {
            Array<T> arr(n);
            for (int i = 0; i < n; i++) arr[i] = Uniform(a, b);
            return arr;
        }
        VectorXd Random_VectorXd(int n, real a = 0.0, real b = 1.0);
        VectorXi Random_VectorXi(int n, int a = 0, int b = 1);
        real Random_Sign(void);

        template<class T>
        void Sparse_Diagonal_Dominant_Matrix(int cols, int rows, Eigen::SparseMatrix<T, Eigen::RowMajor, int>& mat) {
            std::vector<Eigen::Triplet<T>> tripletList;
            tripletList.reserve(3 * rows); //three entries per row
            for (int i = 0; i < rows; i++)
            {
                real sum = 0;
                for (int j = std::max(i - 1, 0); j <= std::min(i + 1, cols - 1); j++) {
                    if (i != j) {
                        real random = Random::Random();
                        sum += abs(random);
                        tripletList.push_back(Eigen::Triplet<T>(i, j, random));
                    }
                }

                sum += abs(Random::Random()); //random offset
                tripletList.push_back(Eigen::Triplet<T>(i, i, 100*sum*Random_Sign()));
            }

            mat.resize(rows, cols);
            mat.setFromTriplets(tripletList.begin(), tripletList.end());
        }
	}

}