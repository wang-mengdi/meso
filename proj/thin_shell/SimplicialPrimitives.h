#pragma once
#include "Common.h"
#include "AuxFunc.h"
#include "GeometryPrimitives.h"

using namespace Meso;
//////////////////////////////////////////////////////////////////////////
////Segment
//////////////////////////////////////////////////////////////////////////

template<int d> class Segment{};

////segment always has a positive phi
template<> class Segment<2>
{Typedef_VectorD(2);
public:
	VectorD p0;
	VectorD p1;

	Segment() :p0(), p1() {}
	Segment(const VectorD& _p0,const VectorD& _p1):p0(_p0),p1(_p1){}

	virtual bool Inside(const VectorD& pos) const {return false;}	////always outside for a 2D segment
	virtual real Phi(const VectorD& pos) const {return Phi(p0,p1,pos);}
	virtual VectorD Normal(const VectorD& pos) const {return Normal(p0,p1,pos);}

	static real Phi(const VectorD& p0,const VectorD& p1,const VectorD& pos)
	{return Closest_Point_Vector(p0,p1,pos).norm();}

	static VectorD Normal(const VectorD& p0,const VectorD& p1,const VectorD& pos)
	{return Closest_Point_Vector(p0,p1,pos).normalized();}

	static VectorD Closest_Point_Vector(const VectorD& p0,const VectorD& p1,const VectorD& pos)
	{real alpha;return Closest_Point_Vector_With_Alpha(p0,p1,pos,alpha);}

	static VectorD Closest_Point_Vector_With_Alpha(const VectorD& p0,const VectorD& p1,const VectorD& pos,real& alpha)
	{
		VectorD m=p1-p0;
		VectorD pq=pos-p0;
		alpha=pq.dot(m)/m.squaredNorm();
		if(alpha<(real)0){alpha=(real)0;return (pos-p0);}
		else if(alpha>(real)1){alpha=(real)1;return (pos-p1);}
		else return (pq-alpha*m);			
	}

	////distance between two segments
	static real Distance(const VectorD& p0,const VectorD& p1,const VectorD& q0,const VectorD& q1)
	{
		Vector<real,4> dis;
		dis[0]=Phi(p0,p1,q0);
		dis[1]=Phi(p0,p1,q1);
		dis[2]=Phi(q0,q1,p0);
		dis[3]=Phi(q0,q1,p1);
		return dis.minCoeff();
	}
};

template<> class Segment<3>
{Typedef_VectorD(3);
public:
	VectorD p0;
	VectorD p1;

	Segment(const VectorD& _p0,const VectorD& _p1):p0(_p0),p1(_p1){}

	virtual bool Inside(const VectorD& pos) const {return false;}	////always outside a 3D segment
	virtual real Phi(const VectorD& pos) const {return Phi(p0,p1,pos);}
	virtual VectorD Normal(const VectorD& pos) const {return Normal(p0,p1,pos);}

	static real Phi(const VectorD& p0,const VectorD& p1,const VectorD& pos)
	{VectorD v=Closest_Point_Vector(p0,p1,pos);return v.norm();}

	static VectorD Normal(const VectorD& p0,const VectorD& p1,const VectorD& pos)
	{VectorD v=Closest_Point_Vector(p0,p1,pos);return v.normalized();}

	static VectorD Closest_Point_Vector(const VectorD& p0,const VectorD& p1,const VectorD& pos)
	{
		real alpha;return Closest_Point_Vector_With_Alpha(p0,p1,pos,alpha);
	}

	static VectorD Closest_Point_Vector_With_Alpha(const VectorD& p0,const VectorD& p1,const VectorD& pos,real& alpha)
	{
		VectorD m=p1-p0;
		VectorD pq=pos-p0;
		alpha=pq.dot(m)/m.squaredNorm();
		if(alpha<(real)0){alpha=(real)0;return (pos-p0);}
		else if(alpha>(real)1){alpha=(real)1;return (pos-p1);}
		else return (pq-alpha*m);			
	}
};

//////////////////////////////////////////////////////////////////////////
////Triangle
//////////////////////////////////////////////////////////////////////////

template<int d> class Triangle{};

template<> class Triangle<2>
{Typedef_VectorD(2);
public:
	VectorD p0,p1,p2;

	Triangle(const VectorD& _p0,const VectorD& _p1,const VectorD& _p2):p0(_p0),p1(_p1),p2(_p2){}
	
	virtual bool Inside(const VectorD& pos) const {return Inside(p0,p1,p2,pos);}
	virtual real Phi(const VectorD& pos) const {return Phi(p0,p1,p2,pos);}
	virtual VectorD Normal(const VectorD& pos) const {return Normal(p0,p1,p2,pos);}

	static bool Inside(const VectorD& p0,const VectorD& p1,const VectorD& p2,const VectorD& pos)
	{
		VectorD v0=p1-p0;VectorD q0=pos-p0;real c0=VectorFunc::Cross(v0,q0)[0];
		VectorD v1=p2-p1;VectorD q1=pos-p1;real c1=VectorFunc::Cross(v1, q1)[0];
		VectorD v2=p0-p2;VectorD q2=pos-p2;real c2=VectorFunc::Cross(v2, q2)[0];
		return (c0>(real)0&&c1>(real)0&&c2>(real)0)||(c0<(real)0&&c1<(real)0&&c2<(real)0);
	}

	static real Phi(const VectorD& p0,const VectorD& p1,const VectorD& p2,const VectorD& pos)
	{
		Vector<real,3> phi;
		phi[0]=Segment<2>::Phi(p0,p1,pos);
		phi[1]=Segment<2>::Phi(p1,p2,pos);
		phi[2]=Segment<2>::Phi(p2,p0,pos);
		real sign=Inside(p0,p1,p2,pos)?(real)-1:(real)1;
		return sign*phi[VectorFunc::Abs_Min_Index(phi)];
	}

	static VectorD Normal(const VectorD& p0,const VectorD& p1,const VectorD& p2,const VectorD& pos)
	{
		VectorD v=Closest_Point_Vector(p0,p1,p2,pos);
		real coef=Inside(p0,p1,p2,pos)?(real)-1:(real)1;
		return coef*v.normalized();
	}

	static VectorD Closest_Point_Vector(const VectorD& p0,const VectorD& p1,const VectorD& p2,const VectorD& pos)
	{
		Vector<VectorD,3> v;
		v[0]=Segment<2>::Closest_Point_Vector(p0,p1,pos);
		v[1]=Segment<2>::Closest_Point_Vector(p1,p2,pos);
		v[2]=Segment<2>::Closest_Point_Vector(p2,p0,pos);
		
		Vector<real,3> abs_phi;
		abs_phi[0]=v[0].norm();
		abs_phi[1]=v[1].norm();
		abs_phi[2]=v[2].norm();

		int idx= VectorFunc::Min_Index(abs_phi);
		return v[idx];
	}

	static VectorD Closest_Point_Vector_With_Barycentric(const VectorD& p0,const VectorD& p1,const VectorD& p2,const VectorD& pos,bool& inside,VectorD& barycentric)
	{
		inside=Inside(p0,p1,p2,pos);
		if(inside){barycentric=MeshFunc::Barycentric_Coord(p0,p1,p2,pos);}

		Vector<VectorD,3> v;Vector<real,3> a;
		v[0]=Segment<2>::Closest_Point_Vector_With_Alpha(p0,p1,pos,a[0]);
		v[1]=Segment<2>::Closest_Point_Vector_With_Alpha(p1,p2,pos,a[1]);
		v[2]=Segment<2>::Closest_Point_Vector_With_Alpha(p2,p0,pos,a[2]);
		
		Vector<real,3> abs_phi;
		abs_phi[0]=v[0].norm();
		abs_phi[1]=v[1].norm();
		abs_phi[2]=v[2].norm();

		int idx= VectorFunc::Min_Index(abs_phi);
		
		if(!inside){	////extrapolate barycentric coordinate
			barycentric[idx]=(real)1-a[idx];
			barycentric[(idx+1)%3]=a[idx];}

		return v[idx];			
	}
};

////triangle 3d always has a positive phi
template<> class Triangle<3>
{Typedef_VectorD(3);
public:
	struct ClosestPointFeature
	{
		enum struct Type : uint8_t {
			Undefined = 0,
			Face = 1,
			Edge = 2,
			Vertex = 3,
		};

		uint8_t m_Index : 4;
		Type m_Type		: 4;

	public:
		ClosestPointFeature(void) : m_Type(Type::Undefined), m_Index(0) {}
		ClosestPointFeature(Type const type, uint8_t const index = 0) : m_Type(type), m_Index(index) {}

		static ClosestPointFeature Face() { return ClosestPointFeature(Type::Face); }
		static ClosestPointFeature Edge(uint8_t index) { return ClosestPointFeature(Type::Edge, index); }
		static ClosestPointFeature Vertex(uint8_t index) { return ClosestPointFeature(Type::Vertex, index); }

		static constexpr unsigned int c_EdgeStartPoint[3] = {0, 1, 2};
		static constexpr unsigned int c_EdgeEndPoint[3]	= {1, 2, 0};
	};

	VectorD p0,p1,p2;

	Triangle(): p0(),p1(),p2() {}
	Triangle(const VectorD& _p0,const VectorD& _p1,const VectorD& _p2):p0(_p0),p1(_p1),p2(_p2){}

	virtual bool Inside(const VectorD& pos) const {return false;}		////always outside a 3D triangle
	virtual real Phi(const VectorD& pos) const {return Phi(p0,p1,p2,pos);}
	virtual VectorD Normal(const VectorD& pos) const {return Normal(p0,p1,p2,pos);}

	static bool Inside(const VectorD& p0,const VectorD& p1,const VectorD& p2,const VectorD& pos) {return false;}

	static bool Projection_Inside(const VectorD& p0,const VectorD& p1,const VectorD& p2,const VectorD& pos)
	{
		VectorD e0=p1-p0;VectorD e1=p2-p1;VectorD e2=p0-p2;
		VectorD n=-(e0.cross(e2)).normalized();
		VectorD q=pos-(pos-p0).dot(n)*n;
		bool inside_e0=(e0.cross(q-p0)).dot(n)>0;
		bool inside_e1=(e1.cross(q-p1)).dot(n)>0;
		bool inside_e2=(e2.cross(q-p2)).dot(n)>0;
		return inside_e0&&inside_e1&&inside_e2;
	}

	static real Phi(const VectorD& p0,const VectorD& p1,const VectorD& p2,const VectorD& pos)
	{VectorD v=Closest_Point_Vector(p0,p1,p2,pos);return v.norm();}

	static VectorD Normal(const VectorD& p0,const VectorD& p1,const VectorD& p2,const VectorD& pos)
	{VectorD v=Closest_Point_Vector(p0,p1,p2,pos);return v.normalized();}

	static VectorD Closest_Point_Vector(const VectorD& p0,const VectorD& p1,const VectorD& p2,const VectorD& pos)
	{
		////projected point inside the triangle
		VectorD e0=p1-p0;VectorD e1=p2-p1;VectorD e2=p0-p2;
		VectorD n=(e0.cross(-e2)).normalized();
		VectorD q=pos-(pos-p0).dot(n)*n;
		bool inside_e0=(e0.cross(q-p0)).dot(n)>0;
		bool inside_e1=(e1.cross(q-p1)).dot(n)>0;
		bool inside_e2=(e2.cross(q-p2)).dot(n)>0;
		if(inside_e0&&inside_e1&&inside_e2)return pos-q;

		////projected point outside the triangle
		Vector<VectorD,3> v;
		v[0]=Segment<3>::Closest_Point_Vector(p0,p1,pos);
		v[1]=Segment<3>::Closest_Point_Vector(p1,p2,pos);
		v[2]=Segment<3>::Closest_Point_Vector(p2,p0,pos);
		
		Vector<real,3> abs_phi;
		abs_phi[0]=v[0].norm();
		abs_phi[1]=v[1].norm();
		abs_phi[2]=v[2].norm();

		int idx= VectorFunc::Min_Index(abs_phi);
		return v[idx];
	}

	static VectorD Closest_Point_Vector(const VectorD& p0,const VectorD& p1,const VectorD& p2,const VectorD& pos,ClosestPointFeature& feature)
	{
		////projected point inside the triangle
		VectorD e0=p1-p0;VectorD e1=p2-p1;VectorD e2=p0-p2;
		VectorD n=(e0.cross(-e2)).normalized();
		VectorD q=pos-(pos-p0).dot(n)*n;
		bool inside_e0=(e0.cross(q-p0)).dot(n)>0;
		bool inside_e1=(e1.cross(q-p1)).dot(n)>0;
		bool inside_e2=(e2.cross(q-p2)).dot(n)>0;
		if(inside_e0&&inside_e1&&inside_e2){feature=ClosestPointFeature::Face();return pos-q;}

		////projected point outside the triangle
		Vector<VectorD,3> v;
		Vector<real,3> alphas;
		v[0]=Segment<3>::Closest_Point_Vector_With_Alpha(p0,p1,pos,alphas[0]);
		v[1]=Segment<3>::Closest_Point_Vector_With_Alpha(p1,p2,pos,alphas[1]);
		v[2]=Segment<3>::Closest_Point_Vector_With_Alpha(p2,p0,pos,alphas[2]);

		Vector<real,3> abs_phi;
		abs_phi[0]=v[0].norm();
		abs_phi[1]=v[1].norm();
		abs_phi[2]=v[2].norm();

		int idx= VectorFunc::Min_Index(abs_phi);
		if(alphas[idx]>0&&alphas[idx]<1){feature=ClosestPointFeature::Edge(idx);}
		else if(alphas[idx]==0){feature=ClosestPointFeature::Vertex(ClosestPointFeature::c_EdgeStartPoint[idx]);}
		else{feature=ClosestPointFeature::Vertex(ClosestPointFeature::c_EdgeEndPoint[idx]);}
		return v[idx];
	}

	static VectorD Normal(const VectorD& p0,const VectorD& p1,const VectorD& p2)
	{VectorD e0=p1-p0; VectorD e1=p2-p0;return e0.cross(e1).normalized();}

	static real Area(const VectorD& p0, const VectorD& p1, const VectorD& p2)
	{
		VectorD e0 = p1 - p0; VectorD e1 = p2 - p0;
		return (real)0.5 * e0.cross(e1).norm();
	}

	real Area() const{
		VectorD e0=p1-p0; VectorD e1=p2-p0;
		return (real) 0.5*e0.cross(e1).norm();
	}

	VectorD Normal() const{
		VectorD e0=p1-p0; VectorD e1=p2-p0;
		return e0.cross(e1).normalized();
	}
};