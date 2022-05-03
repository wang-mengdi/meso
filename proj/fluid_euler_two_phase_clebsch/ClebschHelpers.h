//////////////////////////////////////////////////////////////////////////
// Helper functions for this project.
// Author: Author: Anonymous authors
// This file is part of XXX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////

#pragma once

#include <complex>
#include <cmath>
#include "SPX_Common.h"
using namespace std::complex_literals;

//////////////////////////////////////////////////////////////////////////
////wave function helpers
inline Vector4 C2V(const Vector<C,2>& v) { 
	return Vector4(v[0].real(),v[0].imag(),v[1].real(),v[1].imag()); 
} 
inline Vector<C,2> V2C(const Vector4& v) { 
	Vector<C,2> c;c[0]={v[0],v[1]};c[1]={v[2],v[3]}; return c; 
}
inline Vector3 V2A(const Vector4& v) {
	real a=v[0]; real b=v[1]; real c=v[2]; real d=v[3]; 
	real alpha = std::sqrt(a*a+b*b); real beta = std::atan2(b, a); real gamma = std::atan2(d, c); 
	return Vector3(alpha, beta, gamma);
}
inline Vector3 Q2V(const Vector4& v) {
	AngleAxis a=AngleAxis(v);
	Vector3 q=a.axis()*a.angle(); 
	return q;
}
inline Vector4 Q2V(const Quaternion& q) {
	return Vector4(q.w(),q.x(),q.y(),q.z());
}
inline Quaternion V2Q(const Vector4& v) {
	return Quaternion(&v[0]);
}

//// enforced v = h_bar * vel(k in the paper)
template<int d>
inline Vector<C,2> Vel_To_Psi_C(const Vector<real, d>& vel,const Vector<real, d>& pos) {
	Vector<C,2> psi; psi[0]={1.0/sqrt(1.01),0.}; psi[1]={0.1/sqrt(1.01),0.}; 
	real phase=vel.dot(pos); 
	for(int i=0;i<2;i++) { psi[i]*=exp(1i*phase); }
	return psi;
}

template<int d>
inline Vector<C,2> Vel_To_Psi_C(const Vector<real, d>& vel,const Vector<real, d>& pos,const Vector4& init_psi) {
	Vector<C,2> psi = V2C(init_psi);
	//Normalize(psi);
	real phase=vel.dot(pos); 
	for(int i=0;i<2;i++) { psi[i]*=exp(1i*phase); }
	return psi;
}

template<int d>
inline Vector<C, 2> Vel_To_Psi_C(const real& vel, const int& axis, const Vector<real, d>& pos, const Vector4& init_psi) {
	Vector<C, 2> psi = V2C(init_psi);
	//Normalize(psi);
	real phase = vel*pos[axis];
	for (int i=0;i<2;i++) { psi[i]*=exp(1i*phase); }
	return psi;
}

template<int d>
inline Vector4 Vel_To_Psi(const Vector<real, d>& vel,const Vector<real, d>& pos) { 
	return C2V(Vel_To_Psi_C<d>(vel,pos));
}
