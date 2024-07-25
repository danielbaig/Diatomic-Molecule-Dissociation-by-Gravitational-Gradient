#ifndef BLACKHOLE_H
#define BLACKHOLE_H

struct BlackHole
{
    const double mass{ 0. };
    const double a{ 0. }; // angular momentum to mass ratio-ish
    const double l{ 0. }; // gravitomagnetic monopole moment
    const double charge{ 0. }; // electric charge

    const double l2{ l * l };
    const double a2{ a * a };
    const double charge2{ charge * charge };
};

#endif BLACKHOLE_H