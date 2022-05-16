//////////////////////////////////////////////////////////////////////////
// Basic grid data representation
// Copyright (c) (2022-), Bo Zhu, Mengdi Wang, Yitong Deng
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "Common.h"
#include "AuxFunc.h"
#include <functional>

namespace Meso {
	class AttributeBase {
	public:
		virtual ~AttributeBase() {};
		virtual void Resize(const int size)=0;
	};

	template<class T>
	class Attribute : public AttributeBase {
	private:
		T default_value;
		std::shared_ptr<Array<T>> data_ptr = nullptr;
	public:
		virtual ~Attribute() {}
		Attribute(T def_val) : default_value(def_val) {
			data_ptr = std::make_shared<Array<T>>();
		}

		virtual void Resize(const int size) {
			data_ptr->resize(size, default_value);
		}

		virtual Array<T>& Get_Data() const {
			return *data_ptr;
		}
	};

	class Points {
	public:

		int size = 0; // num elements
		std::map<std::string, std::shared_ptr<AttributeBase>> att_map;

		template<class T>
		void Add_Attribute(const std::string name, T default_v) {
			if (att_map.find(name) == att_map.end()) {
				att_map[name] = std::make_shared<Attribute<T>>(default_v);
			}
			else {
				Error("Error: Duplicate variable: {} encountered in Points.h", name);
			}
		}

		template<class T>
		Array<T>& Get_Attribute(const std::string name) {
			if (att_map.find(name) == att_map.end()) {
				Error("Error: Unfound variable: {} in Points.h", name);
			}
			try {
				Attribute<T>& att = dynamic_cast<Attribute<T>&>(*(att_map[name]));
				return att.Get_Data();
			}
			catch (const std::bad_cast& b) {
				Error("Error: Bad type for variable: {} in Points.h", name);
			}
		}

		template<class T>
		T Get_Entry(const std::string name, const int i) {
			Array<T>& att = Get_Attribute<T>(name);
			if (i >= att.size()) {
				Error("Error: Out Of Bounds for index: {} in Points.h", i);
			}
			return att[i];
		}

		int Size(void);
		void Resize(const int _size);
	};
}