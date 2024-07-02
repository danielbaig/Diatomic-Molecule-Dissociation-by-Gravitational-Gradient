// Diatomic Molecule Dissociation.cpp :
// This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <fstream>
#define _USE_MATH_DEFINES
#include "math.h"

#include "Eigen/Dense"



/*
Cebeci H, Ozdemir N, Sentorun S. Motion of the charged test particles in
Kerr-Newman-Taub-NUT spacetime and analytical solutions.
Physical Review D. 2016 May 15;93(10):104031.
*/


constexpr double epsilon_0{ 8.85e-12 };

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

    double t{ 0. };
    double r{ 0. };
    double phi{ 0. };
    double theta{ 0. };

    double mass2{ mass * mass };
    double L_bar{ L / mass };
    double energy_bar{ energy / mass };
    double charge_bar{ charge / mass };

    double r2{ r * r };
};

double sigma(const BlackHole* BH, Particle* p)
{
    double x{ BH->l + BH->a * cos(p->theta) };
    return p->r2 + x * x;
}

double delta(const BlackHole* BH, Particle* p)
{
    return p->r2 - 2. * BH->mass * p->r + BH->a2 - BH->l2 + BH->charge2;
}

double chi(const BlackHole* BH, Particle* p)
{
    const double sintheta{ sin(p->theta) };
    return BH->a * sintheta * sintheta - 2. * BH->l * cos(p->theta);
}


double P_t(const BlackHole* BH, Particle* p, const double sintheta,
    const double chi_p, const double delta_p, const double ral_sqr,
    const bool isPlaner)
{
    /*
    RHS function for the time component.

    Inputs:
        - const BlackHole* BH: Instance of the black hole.
        - Particle* p: Instance of the particle.
        - const double sintheta: Precomputed sin(theta).
        - const double chi_p: Precomputed chi().
        - const double delta_p: Precomputed delta().
        - const double ral_sqr: r^2 + a^2 + l^2.
        - const bool isPlaner: Whether planer orbit condition is met.


    Outputs:
        - double P_t: RHS function.
    */
    const double x{ p->energy_bar * ral_sqr
        - BH->a * p->L_bar - p->charge_bar * BH->charge * p->r };

    return chi_p * !isPlaner * (p->L_bar - p->energy_bar * chi_p) / (sintheta * sintheta) + x * ral_sqr / delta_p;
}



double P_phi(const BlackHole* BH, Particle* p, const double sintheta,
    const double chi_p, const double delta_p, const double ral_sqr,
    const bool isPlaner)
{
   /*
   RHS function for the phi component.

   Inputs:
       - const BlackHole* BH: Instance of the black hole.
       - Particle* p: Instance of the particle.
       - const double sintheta: Precomputed sin(theta).
       - const double chi_p: Precomputed chi().
       - const double delta_p: Precomputed delta().
       - const double ral_sqr: r^2 + a^2 + l^2.
       - const bool isPlaner: Whether planer orbit condition is met.

   Outputs:
       - double P_phi: RHS function.
   */

    const double x{!isPlaner*(p->L_bar - p->energy_bar*chi_p) / (sintheta*sintheta)};
    return x + (ral_sqr 
        - p->L_bar * BH->a - p->charge_bar * BH->charge * p->r) * BH->a / delta_p;
}




double f_r(const BlackHole* BH, Particle* p, const double sigma_p, const double delta_p,
    const double chi_p, const double ral_sqr, const double dlambda)
{
    /*
    All zeroth and first order differential contributions to the radial motion.

    Inputs:
       - const BlackHole* BH: Instance of the black hole.
       - Particle* p: Instance of the particle.
       - const double delta_p: Precomputed delta().
       - const double chi_p: Precomputed chi().
       - const double delta_p: Precomputed delta().
       - const double ral_sqr: r^2 + a^2 + l^2.
       - const double dlambda: Affine parameter step.

   Outputs:
       - double out: RHS function.
    */
    double t_dot{ (p->t - p->t_prev) / dlambda };
    double r_dot{ (p->r - p->r_prev) / dlambda };
    double phi_dot{ (p->phi - p->phi_prev) / dlambda };
    double theta_dot{ (p->theta - p->theta_prev) / dlambda };

    return 2. * r_dot * r_dot / sigma_p * (p->r * sigma_p / delta_p - p->r - BH->mass * sigma_p / delta_p)
        + 2. * BH->a * r_dot * theta_dot / sigma_p * (BH->l + BH->a * cos(p->theta)) * sin(p->theta)
        + 1. / (sigma_p * sigma_p) * (t_dot - chi_p * phi_dot) * (t_dot - chi_p * phi_dot)
        * (p->r * delta_p / sigma_p - p->r + BH->mass)
        + p->r * delta_p / sigma_p * (r_dot * r_dot / delta_p + theta_dot * theta_dot)
        + r_dot * r_dot / delta_p * (BH->mass - p->r)
        - p->r * sin(p->theta) * sin(p->theta) * delta_p / (sigma_p * sigma_p * sigma_p)
        * (BH->a * t_dot - ral_sqr * phi_dot) * (BH->a * t_dot - ral_sqr * phi_dot)
        + 2. * p->r * delta_p * sin(p->theta) * sin(p->theta) / (sigma_p * sigma_p)
        * (ral_sqr * phi_dot * phi_dot - BH->a * t_dot * phi_dot);

}

double f_theta(const BlackHole* BH, Particle* p, const double sigma_p, const double delta_p,
    const double chi_p, const double ral_sqr, const double dlambda)
{
    /*
    All zeroth and first order differential contributions to the polar motion.

     
    Inputs:
       - const BlackHole* BH: Instance of the black hole.
       - Particle* p: Instance of the particle.
       - const double delta_p: Precomputed delta().
       - const double chi_p: Precomputed chi().
       - const double delta_p: Precomputed delta().
       - const double ral_sqr: r^2 + a^2 + l^2.
       - const double dlambda: Affine parameter step.

   Outputs:
       - double out: RHS function.
    */

    double t_dot{ (p->t - p->t_prev) / dlambda };
    double r_dot{ (p->r - p->r_prev) / dlambda };
    double phi_dot{ (p->phi - p->phi_prev) / dlambda };
    double theta_dot{ (p->theta - p->theta_prev) / dlambda };

    return -2. * p->r * r_dot * theta_dot / sigma_p
        + 2. * BH->a / sigma_p * (BH->l + BH->a * cos(p->theta)) 
        * theta_dot * theta_dot * sin(p->theta)
        - delta_p / (sigma_p * sigma_p * sigma_p) * BH->a * sin(p->theta) 
        * (BH->l + BH->a * cos(p->theta))
        * (t_dot - chi_p * phi_dot) * (t_dot - chi_p * phi_dot)
        + delta_p / (sigma_p * sigma_p) * (t_dot - chi_p * phi_dot)
        * phi_dot * (BH->a * sin(2. * p->theta) + 2. * BH->l * sin(p->theta))
        - BH->a / sigma_p * sin(p->theta) * (BH->l + BH->a * cos(p->theta))
        * (r_dot * r_dot / delta_p + theta_dot * theta_dot)
        + (BH->a * t_dot - ral_sqr * phi_dot) * (BH->a * t_dot - ral_sqr * phi_dot)
        * (sin(2. * p->theta) / (2. * sigma_p * sigma_p)
            + BH->a * sin(p->theta) * sin(p->theta) * sin(p->theta)
            / (sigma_p * sigma_p * sigma_p) * (BH->l + BH->a * cos(p->theta)));
}


void eulerMove(const BlackHole* BH, Particle* p,
    const bool isPlaner, const double dlambda)
{
    /*
    Move forwards one step by the Euler method using the forward difference.

    Inputs:
       - const BlackHole* BH: Instance of the black hole.
       - Particle* p: Instance of the particle.
       - const bool isPlaner: Whether planer orbit condition is met.
       - const double dlambda: Simulation affine parameter step.
    */


    const double sintheta{ sin(p->theta) };
    const double ral_sqr{ p->r2 + BH->a2 + BH->l2 };
    const double sigma_p{ sigma(BH, p) };
    const double chi_p{ chi(BH, p) };
    const double delta_p{ delta(BH, p) };

    double dt{ dlambda * P_t(BH, p, sintheta, chi_p, delta_p, ral_sqr, isPlaner) };
    double dphi{ dlambda * P_phi(BH, p, sintheta, chi_p, delta_p, ral_sqr, isPlaner) };

    double next_r{ 2 * p->r - p->r_prev + dlambda * dlambda * f_r(BH, p, sigma_p, delta_p,
        chi_p,  ral_sqr, dlambda) };
    double next_theta{ 2 * p->theta - p->theta_prev + dlambda * dlambda * f_theta(BH,
        p, sigma_p, delta_p,
        chi_p,  ral_sqr, dlambda) };

    // Perform update.
    p->t_prev = p->t;
    p->t += dt;

    p->r_prev = p->r;
    p->r = next_r;
    p->r2 = p->r * p->r;

    p->phi_prev = p->phi;
    p->phi += dphi;

    p->theta_prev = p->theta;
    p->theta = next_theta;
}


int main()
{
    // Properties file parameters.
    double moleculeLength{};
    unsigned int numParticles{};

    // TODO: This part will need to be changed for multiple particles at a later date.
    double mass1{};
    double electricCharge1{};
    double energy1{};
    double angularMomentum1{};
    double startTime1{};
    double startRadius1{};
    double startPhi1{};
    double startTheta1{};

    double BH_mass{};
    double BH_a{};
    double BH_l{};
    double BH_charge{};

    bool isPlaner{ false };


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

        if (line == "# moleculeLength [double]")
        {
            propertiesFile >> line;
            ++lineCount;
            moleculeLength = std::stod(line);
        }
        else if (line == "# numParticles [unsigned int]")
        {
            propertiesFile >> line;
            ++lineCount;
            numParticles = static_cast<unsigned int>(std::stoi(line));
        }

        // TODO: This part will need to be changed for multiple particles at a later date.
        else if (line == "# particle masses (delimit by commas) [double]")
        {
            propertiesFile >> line;
            ++lineCount;
            mass1 = std::stod(line);
        }
        else if (line == "# particle electric charges (delimit by commas) [double]")
        {
            propertiesFile >> line;
            ++lineCount;
            electricCharge1 = std::stod(line);
        }
        else if (line == "# particle energies (delimit by commas) [double]")
        {
            propertiesFile >> line;
            ++lineCount;
            if (line == "planer")
            {
                isPlaner = true;
                continue;
            }
            energy1 = std::stod(line);
        }
        else if (line == "# particle angular momenta (delimit by commas) [double]")
        {
            propertiesFile >> line;
            ++lineCount;
            angularMomentum1 = std::stod(line);
        }
        else if (line == "# particle start times (delimit by commas) [double]")
        {
            propertiesFile >> line;
            ++lineCount;
            startTime1 = std::stod(line);
        }
        else if (line == "# particle start radii (delimit by commas) [double]")
        {
            propertiesFile >> line;
            ++lineCount;
            startRadius1 = std::stod(line);
        }
        else if (line == "# particle start phis (delimit by commas) [double]")
        {
            propertiesFile >> line;
            ++lineCount;
            startPhi1 = std::stod(line);
        }
        else if (line == "# particle start thetas (delimit by commas) [double]")
        {
            propertiesFile >> line;
            ++lineCount;
            startTheta1 = std::stod(line);
        }

        // Black hole properties
        else if (line == "# Black hole mass [double]")
        {
            propertiesFile >> line;
            ++lineCount;
            BH_mass = std::stod(line);
        }
        else if (line == "# Black hole angular momentum to mass ratio-ish (a) [double]")
        {
            propertiesFile >> line;
            ++lineCount;
            BH_a = std::stod(line);
        }
        else if (line == "# Black hole gravitomagnetic monopole moment (l) [double]")
        {
            propertiesFile >> line;
            ++lineCount;
            BH_l = std::stod(line);
        }
        else if (line == "# Black hole electric charge (Q) [double]")
        {
            propertiesFile >> line;
            ++lineCount;
            BH_charge = std::stod(line);
        }

        else if (line == "")
        {
            continue;
        }
        else
        {
            std::cout << "Unrecognised line in properties file: " << line << '\n';
            throw;
        }


    }



    if (BH_a * BH_a + BH_charge * BH_charge > BH_mass * BH_mass)
    {
        std::cout << "Warning: No horizon present.\n";
    }



    const BlackHole BH{ BH_mass, BH_a, BH_l, BH_charge };
    Particle particle1{ mass1, electricCharge1, angularMomentum1, energy1, startTime1,
        startRadius1, startPhi1, startTheta1, startTime1,
        startRadius1, startPhi1, startTheta1 };




    constexpr double dlambda{ 5e-16 };
    double r_dot0{ 1e+5 };
    double theta_dot0{ 1e+10 };
    // First step must be done manually to account for initial conditions.
    const double sintheta{ sin(particle1.theta) };
    const double ral_sqr{ particle1.r_prev*particle1.r_prev + BH.a2 + BH.l2 };
    const double sigma_p{ sigma(&BH, &particle1) };
    const double chi_p{ chi(&BH, &particle1) };
    const double delta_p{ delta(&BH, &particle1) };

    particle1.t = particle1.t_prev + dlambda * P_t(&BH, &particle1, sintheta,
        chi_p, delta_p, ral_sqr, isPlaner);
    particle1.phi = particle1.phi_prev +  dlambda * P_phi(&BH, &particle1, sintheta,
        chi_p, delta_p, ral_sqr, isPlaner);

    particle1.r = particle1.r_prev + dlambda * r_dot0 + dlambda * dlambda * f_r(&BH,
        &particle1, sigma_p, delta_p,
        chi_p, ral_sqr, dlambda);
    particle1.theta = particle1.theta_prev + dlambda * theta_dot0 + dlambda * dlambda 
        * f_theta(&BH, &particle1, sigma_p, delta_p,
        chi_p, ral_sqr, dlambda);



    std::cout << particle1.r << ',' << particle1.phi << '\n';

    std::cout << "energy: " << energy1 << '\n';
    const double r_Q2{ BH.charge * BH.charge / (4. * M_PI * epsilon_0) };
    const double r_s{ 2. * BH.mass };

    // File setup.
    FILE* coords_TXT = NULL;
    coords_TXT = fopen("../data/coords.txt", "w");

    // Loop variables.
    constexpr unsigned int maxStep{ static_cast<unsigned int>(1e+5) };
    const unsigned int period{ maxStep / 10 };
    unsigned int step{ 0 };


    while (step < maxStep)
    {
        if ((step+1)%period==0)
        {
            std::cout << 100* static_cast<double>(step+1) / static_cast<double>(maxStep) << "%\n";
        }
        // Enter coordinate data into file.
        fprintf(coords_TXT, "%lf,%lf,%lf,%lf\n", particle1.t,
            particle1.r, particle1.phi, particle1.theta);
        

        eulerMove(&BH, &particle1, false, dlambda);


        if (particle1.r < r_s / 2. + sqrt(r_s * r_s / 4. - BH.a2 - r_Q2))
        {
            std::cout << "Entered outer horizon.\n";
            printf("%lf,%lf,%lf,%lf\n", particle1.t,
                particle1.r, particle1.phi, particle1.theta);
            break;
        }



        ++step;
    }
    
    fflush(coords_TXT);



    return 0;
}
