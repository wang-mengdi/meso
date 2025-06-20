#pragma once

#include "Common.h"
#include "nanoflann.hpp"
#include <iostream>
#include <memory>
#include <array>
#include <vector>
#include <functional>

namespace Meso {

    // ================= ArraySlice =====================
    template<class T>
    class ArraySlice {
    public:
        size_t beg_idx = 0, end_idx = 0;
        std::shared_ptr<Array<T> > arr_ptr;
        ArraySlice() {}
        ArraySlice(size_t _beg, size_t _end, Array<T>& arr)
            : beg_idx(_beg), end_idx(_end) {
            arr_ptr = std::make_shared<Array<T> >(arr);
        }
        size_t size(void) const { return end_idx - beg_idx; }
        T& operator[](size_t idx) { return (*arr_ptr)[beg_idx + idx]; }
        const T& operator[](size_t idx) const { return (*arr_ptr)[beg_idx + idx]; }
    };

    // ================ NeighborSearcher =================
    template<int d>
    class NeighborSearcher {
        Typedef_VectorD(d);
        using FilterFunc = std::function<bool(const int)>;

    public:
        Array<Array<int> > search_results;

        NeighborSearcher() { search_results.clear(); }

        virtual void Build_Data(Array<VectorD>& points) = 0;
        virtual std::shared_ptr<Array<VectorD> > Points_Ptr(void) = 0;

        virtual size_t Find_Nbs(const VectorD& pos, const real& radius,
            Array<int>& results, bool append = false) const = 0;
        virtual int Find_Nearest_Nb(const VectorD& pos) const = 0;
        virtual int Find_K_Nearest_Nbs(const VectorD& pos, int k,
            Array<int>& results) const = 0;

        void Update_Points(Array<VectorD>& points);
        void Update_Points(Array<VectorD>& points, FilterFunc& filter_func);
        const Array<int>& Neighbors(const int& idx);
        size_t Find_Nbs(const VectorD& pos, const real& radius,
            FilterFunc& filter_func, Array<int>& results,
            bool append = false) const;
        Array<int> Find_Nbs(const VectorD& pos, const real& radius) const;
        Array<int> Find_Nbs(const VectorD& pos, const real& radius,
            FilterFunc& filter_func) const;
        void Record_All_Neighbors(const real& radius);

        template<int nbsize>
        int Find_Nbs(const VectorD& pos, const real& radius,
            std::array<int, nbsize>& results) const;
    };

    template<int d>
    void NeighborSearcher<d>::Update_Points(Array<VectorD>& points) {
        search_results.clear();
        this->Build_Data(points);
    }
    template<int d>
    void NeighborSearcher<d>::Update_Points(Array<VectorD>& points,
        FilterFunc& filter_func) {
        Array<VectorD> temp_array;
        temp_array.reserve(points.size());
        for (size_t i = 0; i < points.size(); i++)
            if (filter_func((int)i)) temp_array.push_back(points[i]);
        this->Update_Points(temp_array);
    }
    template<int d>
    const Array<int>& NeighborSearcher<d>::Neighbors(const int& idx) {
        Assert(0 <= idx && idx < search_results.size(),
            "NeighborSearcher error: index overflow");
        return search_results[idx];
    }
    template<int d>
    size_t NeighborSearcher<d>::Find_Nbs(const VectorD& pos, const real& radius,
        FilterFunc& filter_func,
        Array<int>& results, bool append) const {
        if (!append) results.clear();
        Array<int> temp_results;
        this->Find_Nbs(pos, radius, temp_results, false);
        size_t num = 0;
        for (size_t i = 0; i < temp_results.size(); i++) {
            if (filter_func(temp_results[i])) {
                num++;
                results.push_back(temp_results[i]);
            }
        }
        return num;
    }
    template<int d>
    Array<int> NeighborSearcher<d>::Find_Nbs(const VectorD& pos,
        const real& radius) const {
        Array<int> temp_results;
        this->Find_Nbs(pos, radius, temp_results, false);
        return temp_results;
    }
    template<int d>
    Array<int> NeighborSearcher<d>::Find_Nbs(const VectorD& pos,
        const real& radius,
        FilterFunc& filter_func) const {
        Array<int> temp_results;
        this->Find_Nbs(pos, radius, filter_func, temp_results, false);
        return temp_results;
    }
    template<int d>
    void NeighborSearcher<d>::Record_All_Neighbors(const real& radius) {
        search_results.clear();
        std::shared_ptr<Array<VectorD> > points_arr = this->Points_Ptr();
        size_t n = points_arr->size();
        search_results.resize(n);
#pragma omp parallel for
        for (int i = 0; i < n; i++) {
            this->Find_Nbs((*points_arr)[i], radius, search_results[i]);
        }
    }
    template<int d>
    template<int nbsize>
    int NeighborSearcher<d>::Find_Nbs(const VectorD& pos, const real& radius,
        std::array<int, nbsize>& results) const {
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

    // ================= PointSetAdapter =================
    template<int d>
    class PointSetAdapter {
        Typedef_VectorD(d);
    public:
        std::shared_ptr<Array<VectorD> > points_ptr;
        void Initialize(Array<VectorD>& arr);
        size_t kdtree_get_point_count() const;
        real kdtree_get_pt(const size_t idx, const size_t dim) const;
        template <class BBOX>
        bool kdtree_get_bbox(BBOX&) const { return false; }
    };
    template<int d>
    void PointSetAdapter<d>::Initialize(Array<VectorD>& arr) {
        points_ptr = std::make_shared<Array<VectorD> >(arr);
    }
    template<int d>
    size_t PointSetAdapter<d>::kdtree_get_point_count() const {
        if (!points_ptr) return 0;
        return points_ptr->size();
    }
    template<int d>
    real PointSetAdapter<d>::kdtree_get_pt(const size_t idx, const size_t dim) const {
        if (!points_ptr) return 0;
        return (*points_ptr)[idx][dim];
    }

    // ================= NeighborKDTree ===================
    template<int d>
    class NeighborKDTree : public NeighborSearcher<d> {
        using Base = NeighborSearcher<d>;
        Typedef_VectorD(d);

    public:
        using Base::Find_Nbs;
        const int max_leaf = 10;
        PointSetAdapter<d> points;
        using my_kd_tree_t = nanoflann::KDTreeSingleIndexAdaptor<
            nanoflann::L2_Simple_Adaptor<real, PointSetAdapter<d>>,
            PointSetAdapter<d>, d, size_t>;
        my_kd_tree_t index;

        NeighborKDTree()
            : index(d, points, nanoflann::KDTreeSingleIndexAdaptorParams(max_leaf)) {}

        virtual void Build_Data(Array<VectorD>& arr);
        virtual std::shared_ptr<Array<VectorD> > Points_Ptr(void);
        virtual size_t Find_Nbs(const VectorD& pos, const real& radius,
            Array<int>& results, bool append = false) const;
        virtual int Find_Nearest_Nb(const VectorD& pos) const;
        virtual int Find_K_Nearest_Nbs(const VectorD& pos, int k,
            Array<int>& results) const;
    };

    template<int d>
    void NeighborKDTree<d>::Build_Data(Array<VectorD>& arr) {
        points.Initialize(arr);
        index.buildIndex();
    }
    template<int d>
    std::shared_ptr<Array<Vector<real,d>> > NeighborKDTree<d>::Points_Ptr(void) {
        return points.points_ptr;
    }
    template<int d>
    size_t NeighborKDTree<d>::Find_Nbs(const VectorD& pos, const real& radius,
        Array<int>& results, bool append) const {
        nanoflann::SearchParams params; params.sorted = false;
        std::vector<std::pair<size_t, real> > ret_matches;
        size_t nMatches = index.radiusSearch(pos.data(), radius * radius, ret_matches, params);
        if (!append) results.clear();
        results.reserve(results.size() + nMatches);
        for (size_t i = 0; i < ret_matches.size(); i++)
            results.push_back((int)ret_matches[i].first);
        return nMatches;
    }
    template<int d>
    int NeighborKDTree<d>::Find_Nearest_Nb(const VectorD& pos) const {
        std::array<size_t, 1> ret_index;
        std::array<real, 1> out_dist_sqr;
        size_t num_results = index.knnSearch(pos.data(), 1, &ret_index[0], &out_dist_sqr[0]);
        if (num_results == 0) {
            std::cerr << "[Error]NeighborKDTree::Find_Nearest_Nb fails for "
                << pos.transpose() << std::endl;
            exit(1);
        }
        return (int)ret_index[0];
    }
    template<int d>
    int NeighborKDTree<d>::Find_K_Nearest_Nbs(const VectorD& pos, int k,
        Array<int>& results) const {
        size_t* ret_index = new size_t[k];
        real* out_dist_sqr = new real[k];
        size_t num_results = index.knnSearch(pos.data(), k, ret_index, out_dist_sqr);
        results.resize(num_results);
        for (size_t i = 0; i < num_results; i++)
            results[i] = (int)ret_index[i];
        delete[] ret_index;
        delete[] out_dist_sqr;
        return (int)num_results;
    }

    // Optional: If you want to instantiate in one place
    // template class NeighborSearcher<2>;
    // template class NeighborSearcher<3>;
    // template class PointSetAdapter<2>;
    // template class PointSetAdapter<3>;
    // template class NeighborKDTree<2>;
    // template class NeighborKDTree<3>;

}  // namespace Meso
