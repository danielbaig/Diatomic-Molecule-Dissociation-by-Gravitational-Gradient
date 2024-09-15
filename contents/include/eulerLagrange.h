#ifndef EULERLAGRANGE_H
#define EULERLAGRANGE_H

#include <iostream>
#include <tuple>
#define _USE_MATH_DEFINES
#include "math.h"

#include <boost/multiprecision/cpp_dec_float.hpp>
using namespace boost::multiprecision;
using cpp_dec_float_n = number<cpp_dec_float<40>>;

#include "blackHole.h"
#include "particle.h"


struct Precalculated{
	cpp_dec_float_n* temp0{ new cpp_dec_float_n[48]{} };
	cpp_dec_float_n* temp1{ new cpp_dec_float_n[27]{} };
	cpp_dec_float_n* temp2{ new cpp_dec_float_n[431]{} };


	const cpp_dec_float_n sigma_p2{};
	const cpp_dec_float_n sigma_p3{};
	const cpp_dec_float_n sigma_p4{};
	const cpp_dec_float_n delta_p2{};
	const cpp_dec_float_n delta_p3{};
	const cpp_dec_float_n delta_p4{};
	const cpp_dec_float_n chi_p2{};
	const cpp_dec_float_n chi_p3{};
	const cpp_dec_float_n chi_p4{};

	const cpp_dec_float_n tdot{};
	const cpp_dec_float_n rdot{};
	const cpp_dec_float_n phidot{};
	const cpp_dec_float_n thetadot{};

	const cpp_dec_float_n tsdot{};
	const cpp_dec_float_n rsdot{};
	const cpp_dec_float_n phisdot{};
	const cpp_dec_float_n thetasdot{};
	const cpp_dec_float_n tsddot{};
	const cpp_dec_float_n rsddot{};
	const cpp_dec_float_n phisddot{};
	const cpp_dec_float_n thetasddot{};

	const cpp_dec_float_n nx{};
	const cpp_dec_float_n ny{};
	const cpp_dec_float_n nz{};
	const cpp_dec_float_n nmag{};
	const cpp_dec_float_n nmag2{};
	const cpp_dec_float_n nmag5{};
	const cpp_dec_float_n nmag6{};
	const cpp_dec_float_n nmag7{};
	const cpp_dec_float_n nmag8{};
	const cpp_dec_float_n nmag12{};
	const cpp_dec_float_n nmag13{};
	const cpp_dec_float_n nmag14{};

	const cpp_dec_float_n ral_sqr{};
	const cpp_dec_float_n sintheta{};
	const cpp_dec_float_n sintheta2{};
	const cpp_dec_float_n sintheta3{};
	const cpp_dec_float_n sintheta4{};
	const cpp_dec_float_n sintheta5{};
	const cpp_dec_float_n sintheta6{};
	const cpp_dec_float_n sintheta7{};

	const cpp_dec_float_n costheta{};
	const cpp_dec_float_n sin2theta{};
	const cpp_dec_float_n sinphi{};
	const cpp_dec_float_n cosphi{};

	const cpp_dec_float_n sinthetas{};
	const cpp_dec_float_n costhetas{};
	const cpp_dec_float_n sinphis{};
	const cpp_dec_float_n cosphis{};
	const cpp_dec_float_n sinthetasdot{};
	const cpp_dec_float_n costhetasdot{};
	const cpp_dec_float_n sinphisdot{};
	const cpp_dec_float_n cosphisdot{};
	const cpp_dec_float_n sinthetas2{};
	const cpp_dec_float_n costhetas2{};
	const cpp_dec_float_n sinphis2{};
	const cpp_dec_float_n cosphis2{};
	const cpp_dec_float_n rsdot2{};
	const cpp_dec_float_n phisdot2{};

Precalculated(const BlackHole* BH,
		Particle* p, Particle* ps,
		const cpp_dec_float_n sigma_p, const cpp_dec_float_n delta_p,
		const cpp_dec_float_n chi_p, const cpp_dec_float_n dlambda);

~Precalculated();
};

cpp_dec_float_n calc_t_next(const BlackHole* BH,
	Particle* p, Particle* ps,
	const cpp_dec_float_n sigma_p, const cpp_dec_float_n delta_p,
	const cpp_dec_float_n chi_p, const Precalculated* precalculated,
	const cpp_dec_float_n dlambda,
	cpp_dec_float_n* temp=nullptr);

cpp_dec_float_n calc_r_next(const BlackHole* BH,
	Particle* p, Particle* ps,
	const cpp_dec_float_n sigma_p, const cpp_dec_float_n delta_p,
	const cpp_dec_float_n chi_p, const Precalculated* precalculated,
	const cpp_dec_float_n dlambda,
	cpp_dec_float_n* temp=nullptr);

cpp_dec_float_n calc_phi_next(const BlackHole* BH,
	Particle* p, Particle* ps,
	const cpp_dec_float_n sigma_p, const cpp_dec_float_n delta_p,
	const cpp_dec_float_n chi_p, const Precalculated* precalculated,
	const cpp_dec_float_n dlambda,
	cpp_dec_float_n* temp=nullptr);

cpp_dec_float_n calc_theta_next(const BlackHole* BH,
	Particle* p, Particle* ps,
	const cpp_dec_float_n sigma_p, const cpp_dec_float_n delta_p,
	const cpp_dec_float_n chi_p, const Precalculated* precalculated,
	const cpp_dec_float_n dlambda,
	cpp_dec_float_n* temp=nullptr);

std::tuple<cpp_dec_float_n, cpp_dec_float_n,	cpp_dec_float_n, cpp_dec_float_n> eulerMoveMathematica(const BlackHole* BH,
	Particle* p, Particle* ps,
	const cpp_dec_float_n sigma_p, const cpp_dec_float_n delta_p,	const cpp_dec_float_n chi_p, const cpp_dec_float_n dlambda);

#endif EULERLAGRANGE_H
