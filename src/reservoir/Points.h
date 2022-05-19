//////////////////////////////////////////////////////////////////////////
// Basic grid data representation
// Copyright (c) (2022-), Bo Zhu, Mengdi Wang, Yitong Deng
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "Common.h"
#include "AuxFunc.h"
#include <functional>

////macros to define helper functions for manipulating particle attributes
#define Setup_Attribute(a, T, def_val)\
	public: T& a(const int i) {return _##a->Get_Entry(i);}\
	Array<T>& a##Ref(){return _##a->Get_Data();}\
	const T& a(const int i) const {return _##a->Get_Entry(i);}\
	const Array<T>& a##Ref() const {return _##a->Get_Data();}\
	protected: std::shared_ptr<Attribute<T>> _##a = std::make_shared<Attribute<T>>(def_val, #a);\
	AttributeRegistrar<T> _##a##_reg = AttributeRegistrar<T>(#a, _##a, att_map);\
	public:\

namespace Meso {
	
	class AttributeBase {
	public:
		virtual ~AttributeBase() {};
		virtual void Resize(const int size) = 0;
		virtual int Size(void) = 0;
	};

	template<class T>
	class Attribute : public AttributeBase {
	private:
		std::string name;
		T default_value;
		std::shared_ptr<Array<T>> data_ptr = nullptr;
	public:

		virtual ~Attribute() {}
		Attribute(T def_val, std::string _name) : 
			default_value(def_val), name(_name) {
			data_ptr = std::make_shared<Array<T>>();
		}

		virtual void Resize(const int size) {
			data_ptr->resize(size, default_value);
		}

		virtual int Size(void) {
			return data_ptr->size();
		}

		virtual Array<T>& Get_Data() const {
			return *data_ptr;
		}

		virtual T& Get_Entry(const int i) const {
			return Get_Data()[i];
		}
	};

	template<class T>
	class AttributeRegistrar {
	public:
		virtual ~AttributeRegistrar() {}
		AttributeRegistrar(std::string name, std::shared_ptr<Attribute<T>> att_ptr,
							std::map<std::string, std::shared_ptr<AttributeBase>>& att_map) {
			if (att_map.find(name) == att_map.end()) {
				att_map[name] = att_ptr;
			}
			else {
				Error("Error: Duplicate variable: {} encountered in Points.h", name);
			}
		}
	};

	class Points {
	public:
		int size = 0; // num elements
		std::map<std::string, std::shared_ptr<AttributeBase>> att_map;

		int Size(void);
		void Resize(const int _size);
	};
}