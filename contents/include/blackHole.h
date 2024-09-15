#ifndef BLACKHOLE_H
#define BLACKHOLE_H

#include <boost/multiprecision/cpp_dec_float.hpp>
using namespace boost::multiprecision;
using cpp_dec_float_n = number<cpp_dec_float<40>>;

struct BlackHole
{
    const cpp_dec_float_n mass{ 0. };
    const cpp_dec_float_n a{ 0. }; // angular momentum to mass ratio-ish
    const cpp_dec_float_n l{ 0. }; // gravitomagnetic monopole moment
    const cpp_dec_float_n charge{ 0. }; // electric charge

    const cpp_dec_float_n l2{ l * l };
    const cpp_dec_float_n l3{ l * l2 };
    const cpp_dec_float_n l4{ l2 * l2 };
    const cpp_dec_float_n l5{ l2 * l3 };
    const cpp_dec_float_n l6{ l3 * l3 };
    const cpp_dec_float_n l7{ l3 * l4 };
    const cpp_dec_float_n l8{ l4 * l4 };
    const cpp_dec_float_n l9{ l4 * l5 };
    const cpp_dec_float_n a2{ a * a };
    const cpp_dec_float_n a3{ a * a2 };
    const cpp_dec_float_n a4{ a2 * a2 };
    const cpp_dec_float_n a5{ a2 * a3 };
    const cpp_dec_float_n a6{ a3 * a3 };
    const cpp_dec_float_n a7{ a3 * a4 };
    const cpp_dec_float_n a8{ a4 * a4 };
    const cpp_dec_float_n a9{ a4 * a5 };
    const cpp_dec_float_n a10{ a5 * a5 };

    const cpp_dec_float_n charge2{ charge * charge };
};

#endif BLACKHOLE_H