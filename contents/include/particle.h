#ifndef PARTICLE_H
#define PARTICLE_H


struct Particle
{
    const double mass{ 0. };
    const double charge{ 0. }; // electric charge
    const double L{ 0. }; // angular momentum
    const double energy{ 0. };

    double t_prev{ 0. };
    double r_prev{ 0. };
    double phi_prev{ 0. };
    double theta_prev{ 0. };

    double t_prevprev{ t_prev };
    double r_prevprev{ r_prev };
    double phi_prevprev{ phi_prev };
    double theta_prevprev{ theta_prev };

    double t{ t_prev };
    double r{ r_prev };
    double phi{ phi_prev };
    double theta{ theta_prev };

    double mass2{ mass * mass };
    double L_bar{ L / mass };
    double energy_bar{ energy / mass };
    double charge_bar{ charge / mass };

    double r2{ r * r };

};


#endif PARTICLE_H