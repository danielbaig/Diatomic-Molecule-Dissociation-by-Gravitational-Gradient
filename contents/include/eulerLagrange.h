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

struct Precalculated{
	const cpp_dec_float_25 sigma_p2{};
	const cpp_dec_float_25 sigma_p3{};
	const cpp_dec_float_25 sigma_p4{};
	const cpp_dec_float_25 delta_p2{};
	const cpp_dec_float_25 delta_p3{};
	const cpp_dec_float_25 delta_p4{};
	const cpp_dec_float_25 chi_p2{};
	const cpp_dec_float_25 chi_p3{};
	const cpp_dec_float_25 chi_p4{};

	const cpp_dec_float_25 tdot{};
	const cpp_dec_float_25 rdot{};
	const cpp_dec_float_25 phidot{};
	const cpp_dec_float_25 thetadot{};

	const cpp_dec_float_25 tsdot{};
	const cpp_dec_float_25 rsdot{};
	const cpp_dec_float_25 phisdot{};
	const cpp_dec_float_25 thetasdot{};
	const cpp_dec_float_25 tsddot{};
	const cpp_dec_float_25 rsddot{};
	const cpp_dec_float_25 phisddot{};
	const cpp_dec_float_25 thetasddot{};

	const cpp_dec_float_25 nx{};
	const cpp_dec_float_25 ny{};
	const cpp_dec_float_25 nz{};
	const cpp_dec_float_25 nmag{};
	const cpp_dec_float_25 nmag2{};
	const cpp_dec_float_25 nmag5{};
	const cpp_dec_float_25 nmag6{};
	const cpp_dec_float_25 nmag7{};
	const cpp_dec_float_25 nmag8{};
	const cpp_dec_float_25 nmag12{};
	const cpp_dec_float_25 nmag13{};
	const cpp_dec_float_25 nmag14{};

	const cpp_dec_float_25 ral_sqr{};
	const cpp_dec_float_25 sintheta{};
	const cpp_dec_float_25 sintheta2{};
	const cpp_dec_float_25 sintheta3{};
	const cpp_dec_float_25 sintheta4{};
	const cpp_dec_float_25 sintheta5{};
	const cpp_dec_float_25 sintheta6{};
	const cpp_dec_float_25 sintheta7{};

	const cpp_dec_float_25 costheta{};
	const cpp_dec_float_25 sin2theta{};
	const cpp_dec_float_25 sinphi{};
	const cpp_dec_float_25 cosphi{};

	const cpp_dec_float_25 sinthetas{};
	const cpp_dec_float_25 costhetas{};
	const cpp_dec_float_25 sinphis{};
	const cpp_dec_float_25 cosphis{};
	const cpp_dec_float_25 sinthetasdot{};
	const cpp_dec_float_25 costhetasdot{};
	const cpp_dec_float_25 sinphisdot{};
	const cpp_dec_float_25 cosphisdot{};
	const cpp_dec_float_25 sinthetas2{};
	const cpp_dec_float_25 costhetas2{};
	const cpp_dec_float_25 sinphis2{};
	const cpp_dec_float_25 cosphis2{};
	const cpp_dec_float_25 rsdot2{};
	const cpp_dec_float_25 phisdot2{};

Precalculated(const BlackHole* BH,
		Particle* p, Particle* ps,
		const cpp_dec_float_25 sigma_p, const cpp_dec_float_25 delta_p,
		const cpp_dec_float_25 chi_p, const cpp_dec_float_25 dlambda);
};

std::tuple<cpp_dec_float_25, cpp_dec_float_25,	cpp_dec_float_25, cpp_dec_float_25> eulerMoveMathematica(const BlackHole* BH,
	Particle* p, Particle* ps,
	const cpp_dec_float_25 sigma_p, const cpp_dec_float_25 delta_p,	const cpp_dec_float_25 chi_p, const cpp_dec_float_25 dlambda);

#endif EULERLAGRANGE_H
