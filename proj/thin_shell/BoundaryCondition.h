#pragma once
#include "Common.h"
#include "Hashtable.h"	
using namespace Meso;
template<int d> class BoundaryConditionMesh
{Typedef_VectorD(d);
public:
	Hashtable<int,VectorD> psi_D_values;	////fixed points
	Hashtable<int,VectorD> forces;			////force

	BoundaryConditionMesh(){}
	void Set_Psi_D(const int& node,const VectorD& displacement=VectorD::Zero())
	{psi_D_values[node]=displacement;}
	void Set_Force(const int& node,const VectorD& force)
	{forces[node]=force;}	
	bool Is_Psi_D(const int& node) const
	{return psi_D_values.find(node)!=psi_D_values.end();}
	void Clear_Force() 
	{forces.clear();}
	void Clear_Psi_D()
	{psi_D_values.clear();}
};
