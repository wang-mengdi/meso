//////////////////////////////////////////////////////////////////////////
// Hashtable
// Copyright (c) (2022-), Bo Zhu and Fan Feng
// This file is part of Meso, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "Common.h"
#include "AuxFunc.h"
#include <unordered_map>
#include <unordered_set>
#include <functional>

namespace Meso {
	template<class T> struct MesoHash {
		typedef T argument_type; typedef std::size_t result_type;
		result_type operator()(argument_type const& arg) const {
			if constexpr (std::is_same<T, int>::value) {
				result_type const h1(std::hash<int>()(arg)); return h1;
			}
			else if constexpr (std::is_same<T, Vector2i>::value) {
				result_type const h1(std::hash<int>()(arg[0])); result_type const h2(std::hash<int>()(arg[1])); return h1 ^ (h2 << 1);
			}
			else if constexpr (std::is_same<T, Vector3i>::value) {
				result_type const h1(std::hash<int>()(arg[0])); result_type const h2(std::hash<int>()(arg[1]));
				result_type const h3(std::hash<int>()(arg[2])); return h1 ^ (h2 << 1) ^ h3;
			}
			else {
				Error("type is not supported");
				return result_type();
			}
		}
	};

	////Hashtable
	template<class T_KEY, class T> using Hashtable = std::unordered_map<T_KEY, T, MesoHash<T_KEY>>;
	template<class T_KEY, class T> using HashtableMultiValue = std::unordered_multimap<T_KEY, T, MesoHash<T_KEY>>;
	////Hashset
	template<class T_KEY> using Hashset = std::unordered_set<T_KEY,MesoHash<T_KEY>>;

	////Hashtable operations
	template<class T_KEY, class T> bool Has(HashtableMultiValue<T_KEY, T>& hashtable, const T_KEY& key, const T& v)
	{
		auto range = hashtable.equal_range(key);
		for (auto iter = range.first; iter != range.second; iter++) { if (iter->second == v) { return true; } }return false;
	}

	template<class T_KEY, class T> bool Add(HashtableMultiValue<T_KEY, T>& hashtable, const T_KEY& key, const T& v, bool check_duplicate = false)
	{
		bool add = !check_duplicate || !Has(hashtable, key, v); if (add)hashtable.insert(std::make_pair(key, v)); return add;
	}

	template<class T_KEY, class T> void Value_Array(const HashtableMultiValue<T_KEY, T>& hashtable, const T_KEY& key, Array<T>& values)
	{
		auto range = hashtable.equal_range(key); for (auto iter = range.first; iter != range.second; iter++) { values.push_back(iter->second); }
	}

	//////////////////////////////////////////////////////////////////////////
	////sorting related functions

	////Vector index sorting by preserving cw/ccw order
	template<int d> Vector<int, d> Reordered(const Vector<int, d>& e, const int v0)
	{
		if constexpr (d == 2) {
			return v0 == e[0] ? Vector<int, 2>(e[0], e[1]) : Vector<int, 2>(e[1], e[0]);
		}
		else if constexpr (d == 3) {
			return v0 == e[0] ? Vector<int, 3>(e[0], e[1], e[2]) : v0 == e[1] ? Vector<int, 3>(e[1], e[2], e[0]) : Vector<int, 3>(e[2], e[0], e[1]);
		}
		else {
			Error("Dimension not allowed");
		}
	}

	////Vector index sorting for a unique key value
	template<int d> Vector<int, d> Unique_Ordered(const Vector<int, d> e) {
		if constexpr (d == 1) {
			return e;
		}
		else if constexpr (d == 2) {
			return MathFunc::Sorted<2>(e);
		}
		else if constexpr (d == 3) {
			return Reordered<3>(e, e.minCoeff());
		}
		else {
			Error("Dimension not allowed");
		}
	}
}