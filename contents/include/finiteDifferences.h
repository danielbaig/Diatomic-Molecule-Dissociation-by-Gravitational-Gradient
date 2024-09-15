#ifndef FINITEDIFFERENCES_H
#define FINITEDIFFERENCES_H


#define _USE_MATH_DEFINES
#include "math.h"

#include <boost/multiprecision/cpp_dec_float.hpp>
using namespace boost::multiprecision;

using cpp_dec_float_n = number<cpp_dec_float<40>>;


#include "blackHole.h"
#include "particle.h"
#include "eulerLagrange.h"


cpp_dec_float_n sigma(const BlackHole* BH, Particle* p);

cpp_dec_float_n delta(const BlackHole* BH, Particle* p);

cpp_dec_float_n chi(const BlackHole* BH, Particle* p);


void updateCoordinates(Particle* p, const cpp_dec_float_n next_t, const cpp_dec_float_n next_r,
    const cpp_dec_float_n next_phi, const cpp_dec_float_n next_theta);



void eulerMove(const BlackHole* BH, Particle* p1, Particle* p2, const cpp_dec_float_n dlambda);


void setupInitialStep(Particle* p, const cpp_dec_float_n next_t, const cpp_dec_float_n next_r,
    const cpp_dec_float_n next_phi, const cpp_dec_float_n next_theta, const cpp_dec_float_n dlambda,
    const cpp_dec_float_n r_dot0, const cpp_dec_float_n r_sintheta_phi_dot0,
    const cpp_dec_float_n r_theta_dot0);




void applyInitialConditions(const BlackHole* BH, Particle* p1, Particle* p2, cpp_dec_float_n dlambda,
    cpp_dec_float_n r_dot0, cpp_dec_float_n r_sintheta_phi_dot0, cpp_dec_float_n r_theta_dot0);




#endif FINITEDIFFERENCES_H