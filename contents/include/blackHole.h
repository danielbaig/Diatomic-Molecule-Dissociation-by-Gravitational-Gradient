#ifndef BLACKHOLE_H
#define BLACKHOLE_H

#include <boost/multiprecision/cpp_dec_float.hpp>
using namespace boost::multiprecision;
using cpp_dec_float_25 = number<cpp_dec_float<25>>;

struct BlackHole
{
    const cpp_dec_float_25 mass{ 0. };
    const cpp_dec_float_25 a{ 0. }; // angular momentum to mass ratio-ish
    const cpp_dec_float_25 l{ 0. }; // gravitomagnetic monopole moment
    const cpp_dec_float_25 charge{ 0. }; // electric charge

    const cpp_dec_float_25 l2{ l * l };
    const cpp_dec_float_25 l3{ l * l2 };
    const cpp_dec_float_25 l4{ l2 * l2 };
    const cpp_dec_float_25 l5{ l2 * l3 };
    const cpp_dec_float_25 l6{ l3 * l3 };
    const cpp_dec_float_25 l7{ l3 * l4 };
    const cpp_dec_float_25 l8{ l4 * l4 };
    const cpp_dec_float_25 l9{ l4 * l5 };
    const cpp_dec_float_25 a2{ a * a };
    const cpp_dec_float_25 a3{ a * a2 };
    const cpp_dec_float_25 a4{ a2 * a2 };
    const cpp_dec_float_25 a5{ a2 * a3 };
    const cpp_dec_float_25 a6{ a3 * a3 };
    const cpp_dec_float_25 a7{ a3 * a4 };
    const cpp_dec_float_25 a8{ a4 * a4 };
    const cpp_dec_float_25 a9{ a4 * a5 };
    const cpp_dec_float_25 a10{ a5 * a5 };

    const cpp_dec_float_25 charge2{ charge * charge };
};

#endif BLACKHOLE_H