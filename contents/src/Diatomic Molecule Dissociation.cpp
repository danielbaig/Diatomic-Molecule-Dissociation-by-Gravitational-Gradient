// Diatomic Molecule Dissociation.cpp :
// This file contains the 'main' function. Program execution begins and ends there.
//

#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>

#define _USE_MATH_DEFINES
#include "math.h"

#include <boost/multiprecision/cpp_dec_float.hpp>
using namespace boost::multiprecision;
constexpr unsigned int precision{ 40 };

using cpp_dec_float_n = number<cpp_dec_float<precision>>;

#include "blackHole.h"
#include "particle.h"
#include "finiteDifferences.h"





/* 
Myers AL.Natural system of units in general relativity. 2016.
https://www.seas.upenn.edu/~amyers/NaturalUnits.pdf
*/
const cpp_dec_float_n epsilon_0{ 1. / (4. * M_PI) }; // Geometrised-Gaussian units.



std::tuple<cpp_dec_float_n, cpp_dec_float_n, cpp_dec_float_n,
    cpp_dec_float_n, cpp_dec_float_n, cpp_dec_float_n,
    cpp_dec_float_n, cpp_dec_float_n, cpp_dec_float_n,
    cpp_dec_float_n, cpp_dec_float_n, cpp_dec_float_n,
    cpp_dec_float_n, cpp_dec_float_n, cpp_dec_float_n,
    cpp_dec_float_n, cpp_dec_float_n> readPropertiesFile()
{
    /*
    Reads the properties file to determine the parameters of the system.
    (I am aware of the number of variables being returned, this has been
    done to ease readability below so other parameters can be easily
    added in the future.)

    Outputs:
        ~ cpp_dec_float_n: 'All the system's parameters.'
    */

    // Properties file parameters.
    cpp_dec_float_n mass{};
    cpp_dec_float_n epsilon{};
    cpp_dec_float_n sigma0{};
    cpp_dec_float_n electricCharge{};
    cpp_dec_float_n startTime1{};
    cpp_dec_float_n startRadius1{};
    cpp_dec_float_n startPhi1{};
    cpp_dec_float_n startTheta1{};

    cpp_dec_float_n r_dot0{};
    cpp_dec_float_n r_sintheta_phi_dot0{};
    cpp_dec_float_n r_theta_dot0{};
    cpp_dec_float_n dphi{};
    cpp_dec_float_n drCoeff{};

    cpp_dec_float_n BH_mass{};
    cpp_dec_float_n BH_a{};
    cpp_dec_float_n BH_l{};
    cpp_dec_float_n BH_charge{};

    // Conversions.
    /*
    Myers AL.Natural system of units in general relativity. 2016.
    https://www.seas.upenn.edu/~amyers/NaturalUnits.pdf
    */
    const cpp_dec_float_n kg_to_m{ 1. / cpp_dec_float_n(1.3466e+27) };
    // https://physics.nist.gov/cgi-bin/cuu/Value?e
    const cpp_dec_float_n eV_to_m{ cpp_dec_float_n(1.602176634e-19)
        / cpp_dec_float_n(1.2102e+44) };
    const cpp_dec_float_n e_to_1{ cpp_dec_float_n(1.602176634e-19)
        / cpp_dec_float_n(3.2735042501e+16) };
    const cpp_dec_float_n s_to_m{ 1. / cpp_dec_float_n(3.3356e-9) };
    const cpp_dec_float_n sol_to_m{ cpp_dec_float_n(1.98855e+30) * kg_to_m };

    // Open file.
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
            // https://physics.nist.gov/cgi-bin/cuu/Value?Rukg
            mass = kg_to_m * cpp_dec_float_n(1.66053906892e-27) * cpp_dec_float_n(line);
        }
        else if (line == "# particle epsilon [double] [meV]")
        {
            propertiesFile >> line;
            ++lineCount;
            epsilon = eV_to_m * 1e-3 * cpp_dec_float_n(line);
        }
        else if (line == "# particle sigma [double] [nm]")
        {
            propertiesFile >> line;
            ++lineCount;
            sigma0 = 1e-9 * cpp_dec_float_n(line);
        }
        else if (line == "# particle electric charge [double] [e]")
        {
            propertiesFile >> line;
            ++lineCount;
            electricCharge = e_to_1 * cpp_dec_float_n(line);
        }
        else if (line == "# particle start time [double] [s]")
        {
            propertiesFile >> line;
            ++lineCount;
            startTime1 = s_to_m * cpp_dec_float_n(line);
        }
        else if (line == "# particle start radius [double] [Gm]")
        {
            propertiesFile >> line;
            ++lineCount;
            startRadius1 = 1e+9 * cpp_dec_float_n(line);
        }
        else if (line == "# particle start phi [double] [rad]")
        {
            propertiesFile >> line;
            ++lineCount;
            startPhi1 = cpp_dec_float_n(line);
        }
        else if (line == "# particle start theta [double] [rad]")
        {
            propertiesFile >> line;
            ++lineCount;
            startTheta1 = cpp_dec_float_n(line);
        }
        else if (line == "# particle r_dot0 [double] [ms^-1 / c]")
        {
            propertiesFile >> line;
            ++lineCount;
            r_dot0 = cpp_dec_float_n(line);
        }
        else if (line == "# particle r*sintheta*phi_dot0 [double] [ms^-1 / c]")
        {
            propertiesFile >> line;
            ++lineCount;
            r_sintheta_phi_dot0 = cpp_dec_float_n(line);
        }
        else if (line == "# particle r*theta_dot0 [double] [ms^-1 / c]")
        {
            propertiesFile >> line;
            ++lineCount;
            r_theta_dot0 = cpp_dec_float_n(line);
        }
        else if (line == "# particle dphi [double] [rad]")
        {
            propertiesFile >> line;
            ++lineCount;
            dphi = cpp_dec_float_n(line);
        }
        else if (line == "# particle dr/moleculeLength [double]")
        {
            propertiesFile >> line;
            ++lineCount;
            drCoeff = cpp_dec_float_n(line);
        }

        // Black hole properties
        else if (line == "# Black hole mass [double] [sol]")
        {
            propertiesFile >> line;
            ++lineCount;
            BH_mass = sol_to_m * cpp_dec_float_n(line);
        }
        else if (line == "# Black hole Kerr parameter (a) [double]")
        {
            propertiesFile >> line;
            ++lineCount;
            BH_a = cpp_dec_float_n(line);
        }
        else if (line == "# Black hole gravitomagnetic monopole moment (l) [double] [A m]")
        {
            propertiesFile >> line;
            ++lineCount;
            BH_l = cpp_dec_float_n(line);
        }
        else if (line == "# Black hole electric charge (Q) [double] [e]")
        {
            propertiesFile >> line;
            ++lineCount;
            BH_charge = e_to_1 * cpp_dec_float_n(line);
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

    return { mass, epsilon, sigma0, electricCharge, startTime1,
        startRadius1, startPhi1, startTheta1, r_dot0,
        r_sintheta_phi_dot0, r_theta_dot0, dphi, drCoeff,
        BH_mass, BH_a, BH_l, BH_charge };
}






int main()
{
    // Properties file parameters.
    cpp_dec_float_n mass{};
    cpp_dec_float_n epsilon{};
    cpp_dec_float_n sigma0{};
    cpp_dec_float_n electricCharge{};
    cpp_dec_float_n startTime1{};
    cpp_dec_float_n startRadius1{};
    cpp_dec_float_n startPhi1{};
    cpp_dec_float_n startTheta1{};

    cpp_dec_float_n r_dot0{};
    cpp_dec_float_n r_sintheta_phi_dot0{};
    cpp_dec_float_n r_theta_dot0{};
    cpp_dec_float_n dphi{};
    cpp_dec_float_n drCoeff{};

    cpp_dec_float_n BH_mass{};
    cpp_dec_float_n BH_a{};
    cpp_dec_float_n BH_l{};
    cpp_dec_float_n BH_charge{};

    std::tie(mass, epsilon, sigma0, electricCharge, startTime1,
        startRadius1, startPhi1, startTheta1, r_dot0,
        r_sintheta_phi_dot0, r_theta_dot0, dphi, drCoeff,
        BH_mass, BH_a, BH_l, BH_charge) = readPropertiesFile();
    

    cpp_dec_float_n moleculeLength{ 
        cpp_dec_float_n(1.12246204831) * sigma0}; // 2^1/6 * sigma
    std::cout << "Molecule length: " << moleculeLength << " m\n";


    // System physicallity warnings.
    if (BH_a * BH_a + BH_charge * BH_charge > BH_mass * BH_mass)
    {
        std::cout << "Warning: No horizon present.\n";
    }

    const cpp_dec_float_n speed0{ sqrt(
        r_dot0 * r_dot0 + r_theta_dot0 * r_theta_dot0
        + r_sintheta_phi_dot0 * r_sintheta_phi_dot0) };


    if (speed0 > 1)
    {
        std::cout << "Warning: Travelling faster than the speed of light.\n";
    }
    else if (speed0 != 1 && mass == 0)
    {
        std::cout << "Warning: Light-like particle is not moving at light speed.\n";
    }

    std::cout << "speed0: " << std::setprecision(precision) << speed0 << "c\n";



    const BlackHole BH{ BH_mass, BH_a, BH_l, BH_charge };
    Particle particle1{ mass, epsilon, sigma0, -electricCharge, startTime1,
        startRadius1, startPhi1, startTheta1 };

    Particle particle2{ mass, epsilon, sigma0, +electricCharge, startTime1,
        startRadius1 + drCoeff*moleculeLength, startPhi1 + dphi, startTheta1 };


    const cpp_dec_float_n dlambda{ 1e-3 }; // Minkowski: 1e-6 BH: 1e-3
    applyInitialConditions(&BH, &particle1, &particle2, dlambda,
        r_dot0, r_sintheta_phi_dot0, r_theta_dot0);




    const cpp_dec_float_n r_Q2{ BH.charge * BH.charge / (4. * M_PI * epsilon_0) };
    const cpp_dec_float_n r_s{ 2. * BH.mass };
    if (r_s * r_s / 4. - BH.a * BH.a - r_Q2 < 0)
    {
        std::cout << "Error: r_s*r_s / 4. - BH.a*BH.a - r_Q2 < 0 has been violated.\n";
        std::cout << r_s * r_s / 4. - BH.a * BH.a - r_Q2 << ">=0\n";
        throw;
    }

    // File setup.
    std::ofstream coords1_TXT("../data/coords1.txt");
    std::ofstream coords2_TXT("../data/coords2.txt");

    // Set precision for output
    coords1_TXT << std::setprecision(precision) << std::fixed;
    coords2_TXT << std::setprecision(precision) << std::fixed;

    // Loop variables.
    constexpr unsigned int maxStep{
        static_cast<unsigned int>(1e+4) }; // Debug 1e+3
    const unsigned int period{ maxStep / 10 };
    unsigned int step{ 0 };
    cpp_dec_float_n lambda{ 0. };
    cpp_dec_float_n separation{};


    while (step < maxStep)
    {
        if ((step+1)%period==0)
        {
            std::cout << 100* static_cast<double>(step+1) 
                / static_cast<double>(maxStep) << "%\n";
        }

        // Test dissociation.
        // https://math.stackexchange.com/questions/833002/distance-between-two-points-in-spherical-coordinates
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



        // Enter coordinate data into file.
        coords1_TXT << lambda << ',' << particle1.t << ',' << particle1.r << ','
            << particle1.phi << ',' << particle1.theta << '\n';
        coords2_TXT << lambda << ',' << particle2.t << ',' << particle2.r << ','
            << particle2.phi << ',' << particle2.theta << '\n';

        eulerMove(&BH, &particle1, &particle2, dlambda);



        ++step;
        lambda += dlambda;
    }
    

    coords1_TXT.close();
    coords2_TXT.close();



    return 0;
}
