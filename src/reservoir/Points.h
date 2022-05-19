//////////////////////////////////////////////////////////////////////////
// Basic grid data representation
// Copyright (c) (2022-), Bo Zhu, Mengdi Wang, Yitong Deng
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "Common.h"
#include "AuxFunc.h"
#include <functional>
#include <any>

////macros to define helper functions for manipulating particle attributes
#define Register_Attribute_Shortcuts(a, T, def_val)\
	public: T& a(const int i) {return this->template Get_Entry<T>(#a, i);}\
	Array<T>& a##Ref(){return this->template Get_Attribute<T>(#a);}\
	const T& a(const int i) const {return this->template Get_Entry<T>(#a, i);}\
	const Array<T>& a##Ref() const {return this->template Get_Attribute<T>(#a);}\
	protected: std::shared_ptr<Attribute<T>> _##a = std::make_shared<Attribute<T>>(def_val, #a);\
	AttributeRegistrar<T> _##a##_reg = AttributeRegistrar<T>(#a, _##a, att_map);\
	public:\

#define Setup_Attribute(a, T, def_val)\
	Register_Attribute_Shortcuts(a, T, def_val);\
	virtual void Init_Attribute_##a(){this->template Add_Attribute<T>(#a, def_val);}\

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
			Info("Fucking shit! {}", name);
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
	};

	template<class T>
	class AttributeRegistrar {
	public:
		virtual ~AttributeRegistrar() {}
		AttributeRegistrar(std::string name, std::shared_ptr<Attribute<T>> att_ptr,
							std::map<std::string, std::shared_ptr<AttributeBase>> att_map) {
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
		std::map<std::string, Array<std::any>&> data_map;

		template<class T>
		void Add_Attribute(const std::string name, T default_v) {
			if (att_map.find(name) == att_map.end()) {
				att_map[name] = std::make_shared<Attribute<T>>(default_v, name);
			}
			else {
				Error("Error: Duplicate variable: {} encountered in Points.h", name);
			}
		}

		template<class T>
		Array<T>& Get_Attribute(const std::string name) const {
			if (att_map.find(name) == att_map.end()) {
				Error("Error: Unfound variable: {} in Points.h", name);
			}
			try {
				Attribute<T>& att = dynamic_cast<Attribute<T>&>(*(att_map.at(name)));
				return att.Get_Data();
			}
			catch (const std::bad_cast& b) {
				Error("Error: Bad type for variable: {} in Points.h", name);
			}
		}

		template<class T>
		T& Get_Entry(const std::string name, const int i) const {
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