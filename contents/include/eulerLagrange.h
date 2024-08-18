#ifndef EULERLAGRANGE_H
#define EULERLAGRANGE_H

#include <iostream>
#include <tuple>
#define _USE_MATH_DEFINES
#include "math.h"

#include <boost/multiprecision/cpp_dec_float.hpp>
using namespace boost::multiprecision;
using cpp_dec_float_25 = number<cpp_dec_float<25>>;

#include "blackHole.h"
#include "particle.h"


cpp_dec_float_25 calc_t_next(const BlackHole* BH,
	Particle* p, Particle* ps,
	const cpp_dec_float_25 sigma_p, const cpp_dec_float_25 delta_p,
	const cpp_dec_float_25 chi_p, const cpp_dec_float_25 dlambda);

cpp_dec_float_25 calc_r_next(const BlackHole* BH,
	Particle* p, Particle* ps,
	const cpp_dec_float_25 sigma_p, const cpp_dec_float_25 delta_p,
	const cpp_dec_float_25 chi_p, const cpp_dec_float_25 dlambda);

cpp_dec_float_25 calc_phi_next(const BlackHole* BH,
	Particle* p, Particle* ps,
	const cpp_dec_float_25 sigma_p, const cpp_dec_float_25 delta_p,
	const cpp_dec_float_25 chi_p, const cpp_dec_float_25 dlambda);

cpp_dec_float_25 calc_theta_next(const BlackHole* BH,
	Particle* p, Particle* ps,
	const cpp_dec_float_25 sigma_p, const cpp_dec_float_25 delta_p,
	const cpp_dec_float_25 chi_p, const cpp_dec_float_25 dlambda);

std::tuple<cpp_dec_float_25, cpp_dec_float_25,	cpp_dec_float_25, cpp_dec_float_25> eulerMoveMathematica(const BlackHole* BH,
	Particle* p, Particle* ps,
	const cpp_dec_float_25 sigma_p, const cpp_dec_float_25 delta_p,	const cpp_dec_float_25 chi_p, const cpp_dec_float_25 dlambda);

#endif EULERLAGRANGE_H
