#include "finiteDifferences.h"

const cpp_dec_float_n epsilon_0{ 1. / (4. * M_PI) }; // Geometrised-Gaussian units.


cpp_dec_float_n sigma(const BlackHole* BH, Particle* p)
{
    const cpp_dec_float_n x{ BH->l + BH->a * cos(p->theta) };
    return p->r2 + x * x;
}

cpp_dec_float_n delta(const BlackHole* BH, Particle* p)
{
    return p->r2 - 2. * BH->mass * p->r + BH->a2 - BH->l2 + BH->charge2;
}

cpp_dec_float_n chi(const BlackHole* BH, Particle* p)
{
    const cpp_dec_float_n sintheta{ sin(p->theta) };
    return BH->a * sintheta * sintheta - 2. * BH->l * cos(p->theta);
}


void updateCoordinates(Particle* p, const cpp_dec_float_n next_t, const cpp_dec_float_n next_r,
    const cpp_dec_float_n next_phi, const cpp_dec_float_n next_theta)
{
    /*
    Update the coordinates to the new values.

    Inputs:
        - Particle* p: Instance of the particle which coordinates have changed.
        - const cpp_dec_float_n next_t: Next time.
        - const cpp_dec_float_n next_r: Next radius.
        - const cpp_dec_float_n next_phi: Next phi.
        - const cpp_dec_float_n next_theta: Next theta.

    */

    // Perform update.
    p->t_prevprev = p->t_prev;
    p->t_prev = p->t;
    p->t = next_t;

    p->r_prevprev = p->r_prev;
    p->r_prev = p->r;
    p->r = next_r;
    p->r2 = p->r * p->r;
    p->r2 = p->r * p->r;
    p->r3 = p->r * p->r2;
    p->r4 = p->r2 * p->r2;
    p->r5 = p->r2 * p->r3;
    p->r6 = p->r3 * p->r3;
    p->r7 = p->r3 * p->r4;

    p->phi_prevprev = p->phi_prev;
    p->phi_prev = p->phi;
    p->phi = next_phi;

    p->theta_prevprev = p->theta_prev;
    p->theta_prev = p->theta;
    p->theta = next_theta;
}


void eulerMove(const BlackHole* BH, Particle* p1, Particle* p2, const cpp_dec_float_n dlambda)
{
    /*
    Move forwards one step by the Euler method using the forward difference.

    Inputs:
       - const BlackHole* BH: Instance of the black hole.
       - Particle* p1: Instance of the first particle.
       - Particle* p2: Instance of the second particle.
       - const double dlambda: Simulation affine parameter step.
    */


    static const cpp_dec_float_n r_Q2{ BH->charge * BH->charge / (4. * M_PI * epsilon_0) };
    static const cpp_dec_float_n r_s{ 2. * BH->mass };


    const cpp_dec_float_n sigma_p1{ sigma(BH, p1) };
    const cpp_dec_float_n chi_p1{ chi(BH, p1) };
    const cpp_dec_float_n delta_p1{ delta(BH, p1) };


    cpp_dec_float_n next_t1{};
    cpp_dec_float_n next_phi1{};
    cpp_dec_float_n next_r1{};
    cpp_dec_float_n next_theta1{};

    const cpp_dec_float_n sigma_p2{ sigma(BH, p2) };
    const cpp_dec_float_n chi_p2{ chi(BH, p2) };
    const cpp_dec_float_n delta_p2{ delta(BH, p2) };


    cpp_dec_float_n next_t2{};
    cpp_dec_float_n next_phi2{};
    cpp_dec_float_n next_r2{};
    cpp_dec_float_n next_theta2{};


    if (p1->r > r_s / 2. + sqrt(r_s * r_s / 4. - BH->a2 - r_Q2))
    {
        std::tie(next_t1, next_r1, next_phi1, next_theta1) = eulerMoveMathematica(BH,
            p1, p2,
            sigma_p1, delta_p1, chi_p1, dlambda);
    }
    if (p2->r > r_s / 2. + sqrt(r_s * r_s / 4. - BH->a2 - r_Q2))
    {
        std::tie(next_t2, next_r2, next_phi2, next_theta2) = eulerMoveMathematica(BH,
            p2, p1,
            sigma_p2, delta_p2, chi_p2, dlambda);
    }
    if (p1->r > r_s / 2. + sqrt(r_s * r_s / 4. - BH->a2 - r_Q2))
    {
        updateCoordinates(p1, next_t1, next_r1, next_phi1, next_theta1);
    }
    if (p2->r > r_s / 2. + sqrt(r_s * r_s / 4. - BH->a2 - r_Q2))
    {
        updateCoordinates(p2, next_t2, next_r2, next_phi2, next_theta2);
    }
}


void setupInitialStep(Particle* p, const cpp_dec_float_n next_t, const cpp_dec_float_n next_r,
    const cpp_dec_float_n next_phi, const cpp_dec_float_n next_theta, const cpp_dec_float_n dlambda,
    const cpp_dec_float_n r_dot0, const cpp_dec_float_n r_sintheta_phi_dot0,
    const cpp_dec_float_n r_theta_dot0)
{
    /*
    Perform the first step incorporating the initial velocity conditions.

    Inputs:
        - Particle* p: Particle to move.
        - const cpp_dec_float_n next_t: Next time.
        - const cpp_dec_float_n next_r: Next radius.
        - const cpp_dec_float_n next_phi: Next phi.
        - const cpp_dec_float_n next_theta: Next theta.
        - const cpp_dec_float_n dlambda: Lambda step.
        - const cpp_dec_float_n r_dot0: Initial radial velocity.
        - const cpp_dec_float_n r_sintheta_phi_dot0: Initial azimuthal velocity.
        - const cpp_dec_float_n r_theta_dot0: Initial polar velocity.
    */

    p->t_prev = p->t_prevprev + dlambda * 1. + next_t
        - (2. * p->t_prev - p->t_prevprev);
    p->phi_prev = p->phi_prevprev
        + dlambda * r_sintheta_phi_dot0 / (p->r_prev * sin(p->theta_prev)) + next_phi
        - (2. * p->phi_prev - p->phi_prevprev);

    p->r_prev = p->r_prevprev + dlambda * r_dot0 + next_r
        - (2. * p->r_prev - p->r_prevprev);
    p->theta_prev = p->theta_prevprev + dlambda * r_theta_dot0 / p->r_prev
        + next_theta - (2. * p->theta_prev - p->theta_prevprev);



    p->t = p->t_prev + dlambda * 1. + next_t
        - (2. * p->t - p->t_prev);
    p->phi = p->phi_prev + dlambda * r_sintheta_phi_dot0 / (p->r * sin(p->theta)) + next_phi
        - (2. * p->phi - p->phi_prev);

    p->r = p->r_prev + dlambda * r_dot0 + next_r - (2. * p->r - p->r_prev);
    p->theta = p->theta_prev + dlambda * r_theta_dot0 / p->r
        + next_theta - (2. * p->theta - p->theta_prev);


    p->r2 = p->r * p->r;
    p->r3 = p->r * p->r2;
    p->r4 = p->r2 * p->r2;
    p->r5 = p->r2 * p->r3;
    p->r6 = p->r3 * p->r3;
    p->r7 = p->r3 * p->r4;
}




void applyInitialConditions(const BlackHole* BH, Particle* p1, Particle* p2, cpp_dec_float_n dlambda,
    cpp_dec_float_n r_dot0, cpp_dec_float_n r_sintheta_phi_dot0, cpp_dec_float_n r_theta_dot0)
{
    /*
    Perform the update for the initial step.

    Inputs:
       - const BlackHole* BH: Instance of the black hole.
       - Particle* p1: Instance of the first particle.
       - Particle* p2: Instance of the second particle.
       - const double dlambda: Simulation affine parameter step.
       - const cpp_dec_float_n r_dot0: Initial radial velocity.
       - const cpp_dec_float_n r_sintheta_phi_dot0: Initial azimuthal velocity.
       - const cpp_dec_float_n r_theta_dot0: Initial polar velocity.
    */

    // For sigma, chi and delta.
    p1->r2 = p1->r * p1->r;
    p2->r2 = p2->r * p2->r;

    // First step must be done manually to account for initial conditions.
    const cpp_dec_float_n sigma_p1{ sigma(BH, p1) };
    const cpp_dec_float_n chi_p1{ chi(BH, p1) };
    const cpp_dec_float_n delta_p1{ delta(BH, p1) };

    const cpp_dec_float_n sigma_p2{ sigma(BH, p2) };
    const cpp_dec_float_n chi_p2{ chi(BH, p2) };
    const cpp_dec_float_n delta_p2{ delta(BH, p2) };



    cpp_dec_float_n next_t1{ p1->t };
    cpp_dec_float_n next_phi1{ p1->phi };
    cpp_dec_float_n next_r1{ p1->r };
    cpp_dec_float_n next_theta1{ p1->theta };

    cpp_dec_float_n next_t2{ p2->t };
    cpp_dec_float_n next_phi2{ p2->phi };
    cpp_dec_float_n next_r2{ p2->r };
    cpp_dec_float_n next_theta2{ p2->theta };

    /*
    std::tie(next_t1, next_r1, next_phi1, next_theta1) = eulerMoveMathematica(BH,
        p1, p2,
        sigma_p1, delta_p1, chi_p1, dlambda);

    std::tie(next_t2, next_r2, next_phi2, next_theta2) = eulerMoveMathematica(BH,
        p2, p1,
        sigma_p2, delta_p2, chi_p2, dlambda);
    */

    setupInitialStep(p1, next_t1, next_r1, next_phi1, next_theta1,
        dlambda, r_dot0, r_sintheta_phi_dot0, r_theta_dot0);

    setupInitialStep(p2, next_t2, next_r2, next_phi2, next_theta2,
        dlambda, r_dot0, r_sintheta_phi_dot0, r_theta_dot0);

}

