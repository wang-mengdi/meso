//////////////////////////////////////////////////////////////////////////
// Common header
// Copyright (c) (2018-), Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once

#include "Common.h"
#include "SpatialHashing.h"
#include "nanoflann.hpp"
#include <iostream>

namespace Meso {
	//It performs like an Array<T>, but with only O(1) space
	//Its size is fixed
	template<class T>
	class ArraySlice {
	public:
		size_t beg_idx = 0, end_idx = 0;//[beg_idx, end_idx), similar to vector::begin() and vetor::end()
		std::shared_ptr<Array<T> > arr_ptr;
		ArraySlice() {}
		ArraySlice(size_t _beg, size_t _end, Array<T>& arr) : beg_idx(_beg), end_idx(_end) { arr_ptr = std::make_shared<Array<T> >(arr); }
		const
			size_t size(void)const { return end_idx - beg_idx; }
		T& operator [] (size_t idx) { return (*arr_ptr)[beg_idx + idx]; }
		const T& operator [] (size_t idx)const { return (*arr_ptr)[beg_idx + idx]; }
	};


	//// Usage:
	//// Initialize (only once) -> Update_Points(frame_1_points) -> queries -> Update_Points(frame_2_points) -> ...
	template<int d>
	class NeighborSearcher {
		Typedef_VectorD(d);
		using FilterFunc = std::function<bool(const int)>;
	public:
		Array<Array<int> > search_results;
	public:
		// Create and initialize with NO points
		NeighborSearcher() { search_results.clear(); }

		//// A Child-Class must implement the following functions:
		// Re-build data structure with points
		virtual void Build_Data(Array<VectorD>& points) = 0;
		virtual std::shared_ptr<Array<VectorD> > Points_Ptr(void) = 0;

		// NOTE: all query functions except for Record_Neighbors should be thread safe. If not, please refer to Mengdi.

		// Store results in Array<int>&results. Clear former results when append==false
		virtual size_t Find_Nbs(const VectorD& pos, const real& radius, Array<int>& results, bool append = false)const = 0;
		virtual int Find_Nearest_Nb(const VectorD& pos)const = 0;
		virtual int Find_K_Nearest_Nbs(const VectorD& pos, int k, Array<int>& results)const = 0;


	public:

		//// Other interfaces:

		// Clear and load points, used for every time step when points are updated
		void Update_Points(Array<VectorD>& points);
		void Update_Points(Array<VectorD>& points, FilterFunc& filter_func);

		const Array<int>& Neighbors(const int& idx) {
			Assert(0 <= idx && idx < search_results.size(), "NeighborSearcher error: index overflow");
			return search_results[idx];
		}
		//Search points in radius and with filter_func()==true
		size_t Find_Nbs(const VectorD& pos, const real& radius, FilterFunc& filter_func, Array<int>& results, bool append = false)const;
		//Pass-by-value style, more convenient for use
		Array<int> Find_Nbs(const VectorD& pos, const real& radius)const;
		Array<int> Find_Nbs(const VectorD& pos, const real& radius, FilterFunc& filter_func)const;
		//Save to search_results
		void Record_All_Neighbors(const real& radius);
		//search points, return with std::array style
		template<int nbsize>
		int Find_Nbs(const VectorD& pos, const real& radius, std::array<int, nbsize>& results)const {
			Array<int> temp_arr;
			this->Find_Nbs(pos, radius, temp_arr, false);
			size_t n = temp_arr.size();
			results[0] = 0;
			for (size_t i = 0; i < n; i++) {
				if (results[0] >= nbsize - 1) {
					std::cerr << "Error: [NeighborSearcher] results full\n";
					return -1;
				}
				results[results[0] + 1] = temp_arr[i];
				results[0]++;
			}
			return results[0];
		}
	};

	template<int d>
	class PointSetAdapter {
		//Declare_Eigen_Types(double, d);
		Typedef_VectorD(d);
	public:
		std::shared_ptr<Array<VectorD> > points_ptr;
		void Initialize(Array<VectorD>& arr);
		// Interface required by nanoflann
		size_t kdtree_get_point_count() const;
		real kdtree_get_pt(const size_t idx, const size_t dim) const;
		template <class BBOX>
		bool kdtree_get_bbox(BBOX& /* bb */) const { return false; }
	};

	template<int d>
	class NeighborKDTree : public NeighborSearcher<d> {
		//Declare_Eigen_Types(double, d);
		//using Base = NeighborSearcher<d>;
		Typedef_VectorD(d);
	public:
		const int max_leaf = 10;/* max leaf */
		PointSetAdapter<d> points;
		using my_kd_tree_t = nanoflann::KDTreeSingleIndexAdaptor<
			nanoflann::L2_Simple_Adaptor<real, PointSetAdapter<d> >,//NOTE: It's actually squared L2 norm
			PointSetAdapter<d>,
			d /* dim */,
			size_t /* accessor type */
		>;
		my_kd_tree_t index;
	public:
		NeighborKDTree() :index(d, points, nanoflann::KDTreeSingleIndexAdaptorParams(max_leaf)) { 
			//index.buildIndex(); 
		}
		virtual void Build_Data(Array<VectorD>& arr);
		virtual std::shared_ptr<Array<VectorD> > Points_Ptr(void) { return points.points_ptr; }
		virtual size_t Find_Nbs(const VectorD& pos, const real& radius, Array<int>& results, bool append = false)const;
		virtual int Find_Nearest_Nb(const VectorD& pos)const;
		virtual int Find_K_Nearest_Nbs(const VectorD& pos, int k, Array<int>& results)const;
	};
}
