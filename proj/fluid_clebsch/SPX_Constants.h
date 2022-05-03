//////////////////////////////////////////////////////////////////////////
// Constants declaration
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __SPXConstants_h__
#define __SPXConstants_h__
#include "SPX_Common.h"

////constants
const real pi=(real)3.14159265358979;
const real two_pi=pi*(real)2;
const real half_pi=pi*(real).5;
const real deg_to_rad=pi/(real)180;
const real one_fourth_pi=pi*(real).25;
const real one_third=(real)1/(real)3;
const real two_thirds=(real)2/(real)3;
const real one_sixth=(real)1/(real)6;
const real sqrt_two=(real)sqrt((real)2);
const real sqrt_three=(real)sqrt((real)3);
const real one_over_pi=(real)1/pi;
const real one_over_sqrt_two=(real)1/sqrt_two;
const real one_over_sqrt_three=(real)1/sqrt_three;

namespace PhysicalUnits {
	const real m = 1.0;
	const real cm = 0.01;
	const real mm = 1e-3;
	const real s = 1.0;
	const real ms = 1e-3;
}

namespace PhysicalConstants {
	const real g = 9.8;
	const real rho_water = 1000;
	const real nu_water = 1e-6;//kinematic viscosity, \appprox 1e-6m^2/s
	const real gamma_water = 72.9 * 1e-3;
}

#endif
