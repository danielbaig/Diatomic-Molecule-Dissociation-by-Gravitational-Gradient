// Diatomic Molecule Dissociation.cpp :
// This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <iomanip>
#include <fstream>
#define _USE_MATH_DEFINES
#include "math.h"

//#include "Eigen/Dense"
#include <boost/multiprecision/cpp_dec_float.hpp>
using namespace boost::multiprecision;
using cpp_dec_float_25 = number<cpp_dec_float<25>>;

#include "blackHole.h"
#include "particle.h"
#include "eulerLagrange.h"




/*
Cebeci H, Ozdemir N, Sentorun S. Motion of the charged test particles in
Kerr-Newman-Taub-NUT spacetime and analytical solutions.
Physical Review D. 2016 May 15;93(10):104031.
*/


const cpp_dec_float_25 epsilon_0{ 1. / (4.*M_PI)}; // Geometrised-Gaussian units.
constexpr int precision{25};


cpp_dec_float_25 sigma(const BlackHole* BH, Particle* p)
{
    cpp_dec_float_25 x{ BH->l + BH->a * cos(p->theta) };
    return p->r2 + x * x;
}

cpp_dec_float_25 delta(const BlackHole* BH, Particle* p)
{
    return p->r2 - 2. * BH->mass * p->r + BH->a2 - BH->l2 + BH->charge2;
}

cpp_dec_float_25 chi(const BlackHole* BH, Particle* p)
{
    const cpp_dec_float_25 sintheta{ sin(p->theta) };
    return BH->a * sintheta * sintheta - 2. * BH->l * cos(p->theta);
}


void updateCoordinates(Particle* p, const cpp_dec_float_25 next_t, const cpp_dec_float_25 next_r,
    const cpp_dec_float_25 next_phi, const cpp_dec_float_25 next_theta)
{
    /*
    Update the coordinates to the new values.

    Inputs:
        - Particle* p: Instance of the particle which coordinates have changed.
        - const cpp_dec_float_25 next_t: Next time.
        - const cpp_dec_float_25 next_r: Next radius.
        - const cpp_dec_float_25 next_phi: Next phi.
        - const cpp_dec_float_25 next_theta: Next theta.

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


void eulerMove(const BlackHole* BH, Particle* p1, Particle* p2, const cpp_dec_float_25 dlambda)
{
    /*
    Move forwards one step by the Euler method using the forward difference.

    Inputs:
       - const BlackHole* BH: Instance of the black hole.
       - Particle* p1: Instance of the first particle.
       - Particle* p2: Instance of the second particle.
       - const double dlambda: Simulation affine parameter step.
    */


    static const cpp_dec_float_25 r_Q2{ BH->charge * BH->charge / (4. * M_PI * epsilon_0) };
    static const cpp_dec_float_25 r_s{ 2. * BH->mass };

    const cpp_dec_float_25 sigma_p1{ sigma(BH, p1) };
    const cpp_dec_float_25 chi_p1{ chi(BH, p1) };
    const cpp_dec_float_25 delta_p1{ delta(BH, p1) };

    cpp_dec_float_25 next_t1{};
    cpp_dec_float_25 next_phi1{};
    cpp_dec_float_25 next_r1{};
    cpp_dec_float_25 next_theta1{};

    const cpp_dec_float_25 sigma_p2{ sigma(BH, p2) };
    const cpp_dec_float_25 chi_p2{ chi(BH, p2) };
    const cpp_dec_float_25 delta_p2{ delta(BH, p2) };

    cpp_dec_float_25 next_t2{};
    cpp_dec_float_25 next_phi2{};
    cpp_dec_float_25 next_r2{};
    cpp_dec_float_25 next_theta2{};


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


void setupInitialStep(Particle* p, const cpp_dec_float_25 next_t, const cpp_dec_float_25 next_r,
    const cpp_dec_float_25 next_phi, const cpp_dec_float_25 next_theta, const cpp_dec_float_25 dlambda,
    const cpp_dec_float_25 r_dot0, const cpp_dec_float_25 phi_dot0, const cpp_dec_float_25 theta_dot0)
{
    /*
    Perform the first step incorporating the initial velocity conditions.

    Inputs:
        - Particle* p: Particle to move.
        - const cpp_dec_float_25 next_t: Next time.
        - const cpp_dec_float_25 next_r: Next radius.
        - const cpp_dec_float_25 next_phi: Next phi.
        - const cpp_dec_float_25 next_theta: Next theta.
        - const cpp_dec_float_25 dlambda: Lambda step.
        - const cpp_dec_float_25 r_dot0: Initial radial velocity.
        - const cpp_dec_float_25 phi_dot0: Initial azimuthal velocity.
        - const cpp_dec_float_25 theta_dot0: Initial polar velocity.
    */

    p->t_prev = p->t_prevprev + dlambda * 1. + next_t
        - (2. * p->t_prev - p->t_prevprev);
    p->phi_prev = p->phi_prevprev + dlambda * phi_dot0 + next_phi
        - (2. * p->phi_prev - p->phi_prevprev);

    p->r_prev = p->r_prevprev + dlambda * r_dot0 + next_r
        - (2. * p->r_prev - p->r_prevprev);
    p->theta_prev = p->theta_prevprev + dlambda * theta_dot0
        + next_theta - (2. * p->theta_prev - p->theta_prevprev);



    p->t = p->t_prev + dlambda * 1. + next_t
        - (2. * p->t - p->t_prev);
    p->phi = p->phi_prev + dlambda * phi_dot0 + next_phi
        - (2. * p->phi - p->phi_prev);

    p->r = p->r_prev + dlambda * r_dot0 + next_r - (2. * p->r - p->r_prev);
    p->theta = p->theta_prev + dlambda * theta_dot0
        + next_theta - (2. * p->theta - p->theta_prev);

    p->r2 = p->r * p->r;
    p->r3 = p->r * p->r2;
    p->r4 = p->r2 * p->r2;
    p->r5 = p->r2 * p->r3;
    p->r6 = p->r3 * p->r3;
    p->r7 = p->r3 * p->r4;
}




void applyInitialConditions(const BlackHole* BH, Particle* p1, Particle* p2, cpp_dec_float_25 dlambda,
    cpp_dec_float_25 r_dot0, cpp_dec_float_25 phi_dot0, cpp_dec_float_25 theta_dot0)
{
    /*
    Perform the update for the initial step.

    Inputs:
       - const BlackHole* BH: Instance of the black hole.
       - Particle* p1: Instance of the first particle.
       - Particle* p2: Instance of the second particle.
       - const double dlambda: Simulation affine parameter step.
       - const cpp_dec_float_25 r_dot0: Initial radial velocity.
       - const cpp_dec_float_25 phi_dot0: Initial azimuthal velocity.
       - const cpp_dec_float_25 theta_dot0: Initial polar velocity.
    */

    // First step must be done manually to account for initial conditions.
    const cpp_dec_float_25 sigma_p1{ sigma(BH, p1) };
    const cpp_dec_float_25 chi_p1{ chi(BH, p1) };
    const cpp_dec_float_25 delta_p1{ delta(BH, p1) };

    const cpp_dec_float_25 sigma_p2{ sigma(BH, p2) };
    const cpp_dec_float_25 chi_p2{ chi(BH, p2) };
    const cpp_dec_float_25 delta_p2{ delta(BH, p2) };


    cpp_dec_float_25 next_t1{};
    cpp_dec_float_25 next_phi1{};
    cpp_dec_float_25 next_r1{};
    cpp_dec_float_25 next_theta1{};

    cpp_dec_float_25 next_t2{};
    cpp_dec_float_25 next_phi2{};
    cpp_dec_float_25 next_r2{};
    cpp_dec_float_25 next_theta2{};


    std::tie(next_t1, next_r1, next_phi1, next_theta1) = eulerMoveMathematica(BH,
        p1, p2,
        sigma_p1, delta_p1, chi_p1, dlambda);

    std::tie(next_t2, next_r2, next_phi2, next_theta2) = eulerMoveMathematica(BH,
        p2, p1,
        sigma_p2, delta_p2, chi_p2, dlambda);

    setupInitialStep(p1, next_t1, next_r1, next_phi1, next_theta1,
        dlambda, r_dot0, phi_dot0, theta_dot0);

    setupInitialStep(p2, next_t2, next_r2, next_phi2, next_theta2,
        dlambda, r_dot0, phi_dot0, theta_dot0);
    
}




int main()
{
    // Properties file parameters.
    cpp_dec_float_25 mass{};
    cpp_dec_float_25 epsilon{};
    cpp_dec_float_25 sigma0{};
    cpp_dec_float_25 electricCharge{};
    cpp_dec_float_25 startTime1{};
    cpp_dec_float_25 startRadius1{};
    cpp_dec_float_25 startPhi1{};
    cpp_dec_float_25 startTheta1{};
    cpp_dec_float_25 r_dot0{};
    cpp_dec_float_25 phi_dot0{};
    cpp_dec_float_25 theta_dot0{};
    cpp_dec_float_25 dphi{};
    cpp_dec_float_25 drCoeff{};

    cpp_dec_float_25 BH_mass{};
    cpp_dec_float_25 BH_a{};
    cpp_dec_float_25 BH_l{};
    cpp_dec_float_25 BH_charge{};

    // Conversions.
    const cpp_dec_float_25 kg_to_m{ 1. / 1.3466e+27 };
    const cpp_dec_float_25 eV_to_m{ 1.602176634e-19 / 1.2102e+44 };
    const cpp_dec_float_25 e_to_1{ 1.602176634e-19 / 5.2909e-19 };
    const cpp_dec_float_25 s_to_m{ 1. / 3.3356e-9 };
    const cpp_dec_float_25 sol_to_m{ 1.98855e+30 * kg_to_m };

    std::ifstream propertiesFile{};
    propertiesFile.open("../data/properties.txt");
    if (!propertiesFile.is_open())
    {
        std::cout << "Error: unable to read file: properties.txt\n";
        throw;
    }

    const unsigned int maxFileSize{ 1000 };
    unsigned int lineCount{ 0 };

    std::string line{};
    while (lineCount < maxFileSize && !propertiesFile.eof())
    {
        std::getline(propertiesFile, line);
        ++lineCount;


        if (line == "# particle mass [double] [u]")
        {
            propertiesFile >> line;
            ++lineCount;
            mass = kg_to_m * 1.66053906892e-27 * std::stod(line);
        }
        else if (line == "# particle epsilon [double] [meV]")
        {
            propertiesFile >> line;
            ++lineCount;
            epsilon = eV_to_m * 1e-3 * std::stod(line);
        }
        else if (line == "# particle sigma [double] [nm]")
        {
            propertiesFile >> line;
            ++lineCount;
            sigma0 = 1e-9 * std::stod(line);
        }
        else if (line == "# particle electric charge [double] [e]")
        {
            propertiesFile >> line;
            ++lineCount;
            electricCharge = e_to_1 * std::stod(line);
        }
        else if (line == "# particle start time [double] [s]")
        {
            propertiesFile >> line;
            ++lineCount;
            startTime1 = s_to_m * std::stod(line);
        }
        else if (line == "# particle start radius [double] [Gm]")
        {
            propertiesFile >> line;
            ++lineCount;
            startRadius1 = 1e+9 * std::stod(line);
        }
        else if (line == "# particle start phi [double] [rad]")
        {
            propertiesFile >> line;
            ++lineCount;
            startPhi1 = std::stod(line);
        }
        else if (line == "# particle start theta [double] [rad]")
        {
            propertiesFile >> line;
            ++lineCount;
            startTheta1 = std::stod(line);
        }
        else if (line == "# particle r_dot0 [double] [ms^-1 / c]")
        {
            propertiesFile >> line;
            ++lineCount;
            r_dot0 = std::stod(line);
        }
        else if (line == "# particle phi_dot0 [double] [rad s^-1]")
        {
            propertiesFile >> line;
            ++lineCount;
            phi_dot0 = std::stod(line);
        }
        else if (line == "# particle theta_dot0 [double] [rad s^-1]")
        {
            propertiesFile >> line;
            ++lineCount;
            theta_dot0 = std::stod(line);
        }
        else if (line == "# particle dphi [double] [rad]")
        {
            propertiesFile >> line;
            ++lineCount;
            dphi = std::stod(line);
        }
        else if (line == "# particle dr/moleculeLength [double]")
        {
            propertiesFile >> line;
            ++lineCount;
            drCoeff = std::stod(line);
        }

        // Black hole properties
        else if (line == "# Black hole mass [double] [sol]")
        {
            propertiesFile >> line;
            ++lineCount;
            BH_mass = sol_to_m * std::stod(line);
        }
        else if (line == "# Black hole Kerr parameter (a) [double]")
        {
            propertiesFile >> line;
            ++lineCount;
            BH_a = std::stod(line);
        }
        else if (line == "# Black hole gravitomagnetic monopole moment (l) [double] [A m]")
        {
            propertiesFile >> line;
            ++lineCount;
            BH_l = std::stod(line);
        }
        else if (line == "# Black hole electric charge (Q) [double] [e]")
        {
            propertiesFile >> line;
            ++lineCount;
            BH_charge = e_to_1 * std::stod(line);
        }

        else if (line == "")
        {
            continue;
        }
        else
        {
            std::cout << "Unrecognised line in properties file:\n" << line << '\n';
            throw;
        }


    }
    cpp_dec_float_25 moleculeLength{ 1.12246204831 * sigma0}; // 2^1/6 * sigma
    std::cout << "Molecule length: " << moleculeLength << '\n';


    if (BH_a * BH_a + BH_charge * BH_charge > BH_mass * BH_mass)
    {
        std::cout << "Warning: No horizon present.\n";
    }

    const cpp_dec_float_25 speed0{ r_dot0 * r_dot0 + startRadius1 * startRadius1 * (theta_dot0 * theta_dot0
        + sin(startTheta1) * sin(startTheta1) * phi_dot0 * phi_dot0) };


    if (speed0 > 1)
    {
        std::cout << "Warning: Travelling faster than the speed of light.\n";
    }
    else if (speed0 != 1 && mass == 0)
    {
        std::cout << "Warning: Light-like particle is not moving at light speed.\n";
    }

    std::cout << "speed0: " << speed0 << "c\n";



    const BlackHole BH{ BH_mass, BH_a, BH_l, BH_charge };
    Particle particle1{ mass, epsilon, sigma0, electricCharge, startTime1,
        startRadius1, startPhi1, startTheta1 };

    Particle particle2{ mass, epsilon, sigma0, electricCharge, startTime1,
        startRadius1 + drCoeff*moleculeLength, startPhi1 + dphi, startTheta1 };



    const cpp_dec_float_25 dlambda{ 1e-3 }; // Minkowski: 1e-6 bg: 1e-3
    applyInitialConditions(&BH, &particle1, &particle2, dlambda, r_dot0, phi_dot0, theta_dot0);






    const cpp_dec_float_25 r_Q2{ BH.charge * BH.charge / (4. * M_PI * epsilon_0) };
    const cpp_dec_float_25 r_s{ 2. * BH.mass };
    if (r_s * r_s / 4. - BH.a * BH.a - r_Q2 < 0)
    {
        std::cout << "Error: r_s*r_s / 4. - BH.a*BH.a - r_Q2 < 0 has been violated.\n";
        std::cout << r_s * r_s / 4. - BH.a * BH.a - r_Q2 << ">=0\n";
        throw;
    }

    // File setup.
    /*
    FILE* coords1_TXT = NULL;
    coords1_TXT = fopen("../data/coords1.txt", "w");
    FILE* coords2_TXT = NULL;
    coords2_TXT = fopen("../data/coords2.txt", "w");
    */
    std::ofstream coords1_TXT("../data/coords1.txt");
    std::ofstream coords2_TXT("../data/coords2.txt");

    // Set precision for output
    coords1_TXT << std::setprecision(25) << std::fixed;
    coords2_TXT << std::setprecision(25) << std::fixed;

    // Loop variables.
    constexpr unsigned int maxStep{ static_cast<unsigned int>(1e+4) }; // Debug 1e+3
    const unsigned int period{ maxStep / 10 };
    unsigned int step{ 0 };
    cpp_dec_float_25 lambda{ 0 };
    cpp_dec_float_25 separation{};


    while (step < maxStep)
    {


        if ((step+1)%period==0)
        {
            std::cout << 100* static_cast<double>(step+1) / static_cast<double>(maxStep) << "%\n";
        }
        // Enter coordinate data into file.
        /*
        fprintf(coords1_TXT, "%.25lf,%.25lf,%.25lf,%.25lf,%.25lf\n", lambda, particle1.t,
            particle1.r, particle1.phi, particle1.theta);
        fprintf(coords2_TXT, "%.25lf,%.25lf,%.25lf,%.25lf,%.25lf\n", lambda, particle2.t,
            particle2.r, particle2.phi, particle2.theta);
        */

        if (abs(particle1.phi) > 10 || abs(particle2.phi) > 10)
        {
            std::cout << "Error: phi has become unusually large.\n";
            break;
        }

        separation = sqrt(particle1.r * particle1.r + particle2.r * particle2.r
            - 2. * particle1.r * particle2.r
            * (sin(particle1.theta) * sin(particle2.theta) 
                * cos(particle1.phi - particle2.phi)
                + cos(particle1.theta) * cos(particle2.theta)));

        if (separation > 3. * moleculeLength)
        {
            std::cout << "Possible dissocitation.\n";
            break;
        }

        coords1_TXT << lambda << ',' << particle1.t << ',' << particle1.r << ','
            << particle1.phi << ',' << particle1.theta << '\n';
        coords2_TXT << lambda << ',' << particle2.t << ',' << particle2.r << ','
            << particle2.phi << ',' << particle2.theta << '\n';

        eulerMove(&BH, &particle1, &particle2, dlambda);

        /*

        if (particle1.r < r_s / 2. + sqrt(r_s * r_s / 4. - BH.a2 - r_Q2))
        {
            std::cout << "Entered outer horizon.\n";
            printf("%lf,%lf,%lf,%lf,%lf\n", lambda, particle1.t,
                particle1.r, particle1.phi, particle1.theta);
            std::cout << r_s / 2. + sqrt(r_s * r_s / 4. - BH.a2 - r_Q2) << '\n';
            break;
        }
        else if (particle1.r > 1.1*startRadius1)
        {
            std::cout << "Particle has escaped.\n";
            printf("%lf,%lf,%lf,%lf,%lf\n", lambda, particle1.t,
                particle1.r, particle1.phi, particle1.theta);
            break;
        }
        */


        ++step;
        lambda += dlambda;
    }
    
    // fflush(coords1_TXT);
    // fflush(coords2_TXT);

    coords1_TXT.close();
    coords2_TXT.close();



    return 0;
}
