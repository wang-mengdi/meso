//////////////////////////////////////////////////////////////////////////
// Multigrid method
// Copyright (c) (2022-), Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "LinearMapping.h"

namespace Meso {

	template<class T, int d>
	class VCycleMultigrid :LinearMapping<T> {
		Typedef_VectorD(d);
	public:
		int L;

		Array<LinearMapping<T>*> mappings;
		Array<LinearMapping<T>*> restrictors;
		Array<LinearMapping<T>*> prolongators;
		Array<LinearMapping<T>*> presmoothers;
		Array<LinearMapping<T>*> postsmoothers;

		ArrayDv<ArrayDv<T>> xs;
		ArrayDv<ArrayDv<T>> bs;
		ArrayDv<ArrayDv<T>> rs;
	public:
		virtual int XDof() const {

		}

		virtual int YDof() const {

		}

		//solve some Ax=b
		virtual void Apply(ArrayDv<T>& x, const ArrayDv<T>& b) {

		}

		void Init(void) {

		}

		void V_Cycle(ArrayDv<T>& x0, const ArrayDv<T>& b0) {
			
			for (int i = 0; i < L; i++) {
				presmoothers[i]->
			}
		}
	};

}