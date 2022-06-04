//////////////////////////////////////////////////////////////////////////
// SPH Kernels
// Copyright (c) (2018-),Xiangxin Kong, Mengdi Wang
// Please see simplex/docs/kernels-math-en.md for documentation.
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#include "Kernel.h"

namespace Meso {
	real Kernel::Weight(int d, real r, real h1, KernelType kernel_type)const {
		if (kernel_type == KernelType::NONE) kernel_type = ref_type;
		return kernels[(int)kernel_type]->Weight(d, fabs(r / h1)) / MathFunc::Quick_Pow(h1, d);
	}

	real Kernel::Weight(int d, real r, KernelType kernel_type)const {
		if (kernel_type == KernelType::NONE) kernel_type = ref_type;
		return kernels[(int)kernel_type]->Weight(d, fabs(r / h)) * h_pows_inv[d];
	}

	UnitPOLY6 Kernel::poly6;
	UnitSPIKY Kernel::spiky;
	UnitCUBIC Kernel::cubic;
	UnitQUINTIC Kernel::quintic;
	UnitGAUSSIAN Kernel::gaussian;
}