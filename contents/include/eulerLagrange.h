#ifndef EULERLAGRANGE_H
#define EULERLAGRANGE_H

#include <iostream>
#include <tuple>
#define _USE_MATH_DEFINES
#include "math.h"

#include "blackHole.h"
#include "particle.h"


std::tuple<double, double, double, double> eulerMoveMathematica(const BlackHole* BH,
	Particle* p, Particle* ps,
	const double sigma_p, const double delta_p, const double chi_p, const double dlambda);

#endif EULERLAGRANGE_H
