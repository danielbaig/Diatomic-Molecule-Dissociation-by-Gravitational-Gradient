#ifndef PARTICLE_H
#define PARTICLE_H


#include <boost/multiprecision/cpp_dec_float.hpp>
using namespace boost::multiprecision;
using cpp_dec_float_25 = number<cpp_dec_float<25>>;


struct Particle
{
    const cpp_dec_float_25 mass{ 0. };
    const cpp_dec_float_25 epsilon{ 0. };
    const cpp_dec_float_25 sigma{ 0. };
    const cpp_dec_float_25 charge{ 0. }; // electric charge

    cpp_dec_float_25 t_prev{ 0. };
    cpp_dec_float_25 r_prev{ 0. };
    cpp_dec_float_25 phi_prev{ 0. };
    cpp_dec_float_25 theta_prev{ 0. };

    cpp_dec_float_25 t_prevprev{ t_prev };
    cpp_dec_float_25 r_prevprev{ r_prev };
    cpp_dec_float_25 phi_prevprev{ phi_prev };
    cpp_dec_float_25 theta_prevprev{ theta_prev };

    cpp_dec_float_25 t{ t_prev };
    cpp_dec_float_25 r{ r_prev };
    cpp_dec_float_25 phi{ phi_prev };
    cpp_dec_float_25 theta{ theta_prev };

    cpp_dec_float_25 mass2{ mass * mass };
    cpp_dec_float_25 charge_bar{ charge / mass };

    cpp_dec_float_25 r2{ r * r };
    cpp_dec_float_25 r3{ r * r2 };
    cpp_dec_float_25 r4{ r2 * r2 };
    cpp_dec_float_25 r5{ r2 * r3 };
    cpp_dec_float_25 r6{ r3 * r3 };
    cpp_dec_float_25 r7{ r3 * r4 };

    const cpp_dec_float_25 sigma6{ sigma * sigma * sigma * sigma * sigma * sigma };
    const cpp_dec_float_25 sigma12{ sigma6 * sigma6 };

};


#endif PARTICLE_H