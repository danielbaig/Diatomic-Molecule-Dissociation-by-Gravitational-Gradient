#ifndef PARTICLE_H
#define PARTICLE_H


#include <boost/multiprecision/cpp_dec_float.hpp>
using namespace boost::multiprecision;
using cpp_dec_float_n = number<cpp_dec_float<40>>;


struct Particle
{
    const cpp_dec_float_n mass{ 0. };
    const cpp_dec_float_n epsilon{ 0. };
    const cpp_dec_float_n sigma{ 0. };
    const cpp_dec_float_n charge{ 0. }; // electric charge

    cpp_dec_float_n t_prev{ 0. };
    cpp_dec_float_n r_prev{ 0. };
    cpp_dec_float_n phi_prev{ 0. };
    cpp_dec_float_n theta_prev{ 0. };

    cpp_dec_float_n t_prevprev{ t_prev };
    cpp_dec_float_n r_prevprev{ r_prev };
    cpp_dec_float_n phi_prevprev{ phi_prev };
    cpp_dec_float_n theta_prevprev{ theta_prev };

    cpp_dec_float_n t{ t_prev };
    cpp_dec_float_n r{ r_prev };
    cpp_dec_float_n phi{ phi_prev };
    cpp_dec_float_n theta{ theta_prev };

    cpp_dec_float_n mass2{ mass * mass };
    cpp_dec_float_n charge_bar{ charge / mass };

    cpp_dec_float_n r2{ r * r };
    cpp_dec_float_n r3{ r * r2 };
    cpp_dec_float_n r4{ r2 * r2 };
    cpp_dec_float_n r5{ r2 * r3 };
    cpp_dec_float_n r6{ r3 * r3 };
    cpp_dec_float_n r7{ r3 * r4 };

    const cpp_dec_float_n sigma6{ sigma * sigma * sigma * sigma * sigma * sigma };
    const cpp_dec_float_n sigma12{ sigma6 * sigma6 };

};


#endif PARTICLE_H