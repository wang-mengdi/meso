#pragma once
#include "Common.h"
#include "Constants.h"

using namespace Meso;
namespace ElasticParamImpl 
{
    template<int dim> inline real Lambda(const real youngs,const real poisson){return youngs*poisson/(((real)1+poisson)*((real)1-(real)2*poisson));}
    template<> inline real Lambda<2>(const real youngs,const real poisson){return youngs*poisson/((real)1-pow(poisson,2));}

	template<int dim> inline real Bulk(const real youngs,const real poisson){return (real).5*youngs/((real)1-poisson);}
	template<> inline real Bulk<3>(const real youngs,const real poisson){return CommonConstants::one_third*youngs/((real)1-(real)2*poisson);}

	template<int dim>  inline real E_Hat_Coef(const real poisson){return (real)0;}	////E=E_hat*coef
	template<> inline real E_Hat_Coef<2>(const real poisson){return ((real)1-poisson*poisson);}
	template<> inline real E_Hat_Coef<3>(const real poisson){return ((real)1-2*poisson)*((real)1+poisson);}

	template<int dim> inline real G_Hat_Coef(const real poisson){return (real)0;}	////G=E_hat*coef
	template<> inline real G_Hat_Coef<2>(const real poisson){return (real).5*((real)1-poisson);}	
	template<> inline real G_Hat_Coef<3>(const real poisson){return (real).5*((real)1-(real)2*poisson);}
}

class ElasticParam
{public:
	real youngs_modulus=(real)1.;
	real poisson_ratio=(real).45;
	
	ElasticParam(const real _youngs,const real _poisson):youngs_modulus(_youngs),poisson_ratio(_poisson){}

	////Lame parameters
	template<int dim> real Lambda() const {return Lambda<dim>(youngs_modulus,poisson_ratio);}
	template<int dim> real Mu() const {return Mu<dim>(youngs_modulus,poisson_ratio);}
	
	template <int dim> static real Lambda(const real youngs,const real poisson) {return ElasticParamImpl::Lambda<dim>(youngs, poisson);}
	template<int dim> static real Mu(const real youngs,const real poisson){return (real).5*youngs/((real)1+poisson);}	////shear

	template<int dim> real Bulk() const {return Bulk<dim>(youngs_modulus,poisson_ratio);}
	template<int dim> static real Bulk(const real youngs,const real poisson){return ElasticParamImpl::Bulk<dim>(youngs, poisson);}
	template<int dim> real Shear() const{return Shear<dim>(youngs_modulus,poisson_ratio);}
	template<int dim> static real Shear(const real youngs,const real poisson){return (real).5*youngs/((real)1+poisson);}

	////E_hat
	template<int dim> static real E_Hat_Coef(const real poisson){return ElasticParamImpl::E_Hat_Coef<dim>(poisson);}
	////G_hat
	template<int dim> static real G_Hat_Coef(const real poisson){return ElasticParamImpl::G_Hat_Coef<dim>(poisson);}
};
