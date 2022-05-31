#pragma once
#include "Points.h"

namespace Meso {
	template<int d>
	class DiscreteShellParticles : public Points {
		Typedef_VectorD(d);
	public:

		DiscreteShellParticles() {}

		Setup_Attribute(x, VectorD, VectorD::Zero());
		ArrayPtr<VectorD> xPtr() {return _x->Get_Ptr();}

		Setup_Attribute(v, VectorD, VectorD::Zero());
		Setup_Attribute(f, VectorD, VectorD::Zero());
	};
}