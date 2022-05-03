//////////////////////////////////////////////////////////////////////////
// A class that can point to any Array<T>, to simplify class Particles<d>
// Copyright (c) (2018-), Mengdi Wang
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __ArrayPointer_h__
#define __ArrayPointer_h__
#include <iostream>
#include "SPX_Common.h"
#include "ArrayIO.h"
#include <set>

class ArrayPointerBase {
public:
	virtual ~ArrayPointerBase(){}
	virtual void Resize(const int size) = 0;
	virtual void Reserve(const int size) = 0;
	virtual int Add_Element() = 0;//push back one element. uses push_back instead of resize
	virtual int Add_Elements(const int n) = 0;
	virtual void Copy_Element(const int idx, const std::shared_ptr<ArrayPointerBase> src, const int src_idx) = 0;
	virtual int Delete_Elements(const Array<int>& is_deleted) = 0;//return size after deletion. is_deleted[i]!=0 means a deletion.
	virtual int Delete_Elements(const std::set<int>& to_delete) = 0;//return size after deletion.
	virtual int Read_Binary(const std::string& file_name) = 0;//write size(), then content.
	virtual int Write_Binary(const std::string& file_name) = 0;//pairing with Read_Binary()
	virtual void Read_Binary_Content(std::istream& input, const uint32_t& n) = 0;//do not write size(), just content.
	virtual void Write_Binary_Content(std::ostream& output) = 0;//pairing with Read_Binary_Content()
};

template<class T>
class ArrayPointerDerived :public ArrayPointerBase {
public:
	ArrayPtr<T> array_ptr;
	virtual ~ArrayPointerDerived(){}
	virtual void Resize(const int size) { if (size == 0)array_ptr->clear(); else array_ptr->resize(size, Zero<T>()); }
	virtual void Reserve(const int size) { array_ptr->reserve(size); }
	virtual int Add_Element(void) { array_ptr->push_back(T()); return (int)array_ptr->size() - 1; }
	virtual int Add_Elements(const int n) { array_ptr->resize(array_ptr->size() + n, Zero<T>()); return (int)array_ptr->size() - n; }
	virtual void Copy_Element(const int idx, const std::shared_ptr<ArrayPointerBase> src, const int src_idx) {
		auto derived_src = std::dynamic_pointer_cast<ArrayPointerDerived<T>>(src);
		const ArrayPtr<T> derived_ptr = derived_src->array_ptr;
		(*array_ptr)[idx] = (*derived_ptr)[src_idx];
	}
	virtual int Delete_Elements(const Array<int>& is_deleted) {
		int deleted_size = 0;
		for (int i = 0; i < array_ptr->size(); i++) {
			if (is_deleted[i]) continue;
			(*array_ptr)[deleted_size] = (*array_ptr)[i];
			deleted_size++;
		}
		array_ptr->resize(deleted_size);
		return deleted_size;
	}
	virtual int Delete_Elements(const std::set<int>& to_delete) {
		int deleted_size = 0;
		for (int i = 0; i < array_ptr->size(); i++) {
			if (to_delete.count(i)) continue;
			(*array_ptr)[deleted_size] = (*array_ptr)[i];
			deleted_size++;
		}
		array_ptr->resize(deleted_size);
		return deleted_size;
	}
	virtual int Read_Binary(const std::string& file_name) { return BinaryDataIO::Read_Array(file_name, *array_ptr); }
	virtual int Write_Binary(const std::string& file_name) { return BinaryDataIO::Write_Array(file_name, *array_ptr); }
	virtual void Read_Binary_Content(std::istream& input, const uint32_t& n) { BinaryDataIO::Read_Array_Stream_Content(input, *array_ptr, n); }
	virtual void Write_Binary_Content(std::ostream& output) { BinaryDataIO::Write_Array_Stream_Content(output, *array_ptr); }
	ArrayPointerDerived(ArrayPtr<T> _ptr) :array_ptr(_ptr) { }
};

#endif
