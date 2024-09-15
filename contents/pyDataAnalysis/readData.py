from mpmath import mpmathify, mp

from pathlib import Path
pathtohere = Path.cwd()

mp.dps = 40  # Decimal places of precision

class BlackHole:
    mass = 0.
    a = 0.
    l = 0.
    charge = 0.
    
class Particle:
    mass = 0.
    charge = 0. # electric charge
    energy = 0.
    L = 0. # angular momentum


def readPropertiesFile():
    """
    Read the properties file and save the variables.
    
    Outputs:
        - BH: Instance of the black hole.
        - particle: Instance of the particle.
        - moleculeLength:float: The average molecule length.
    """
    
    BH = BlackHole()
    particle = Particle()
    
    kg_to_m = 1. / mpmathify(1.3466e+27)
    eV_to_m = mpmathify(1.602176634e-19) / mpmathify(1.2102e+44)
    e_to_1 = mpmathify(1.602176634e-19) / mpmathify(5.2909e-19)
    s_to_m = 1. / mpmathify(3.3356e-9)
    sol_to_m = mpmathify(1.98855e+30) * kg_to_m
    
    line = True
    filePath = pathtohere / 'data/properties.txt'
    i = 0
    maxFileSize = 1_000
    # https://stackoverflow.com/questions/3277503/how-to-read-a-file-line-by-line-into-a-list
    with open(filePath, 'r', encoding='UTF-8') as f:
        while line and i < maxFileSize:
            line = f.readline()
            

            if (line[:-1] == "# particle mass [double] [u]"):
                line = f.readline()
                i += 1;
                particle.mass = kg_to_m*mpmathify(1.66053906892e-27)*mpmathify(line)
            elif (line[:-1] == "# particle epsilon [double] [meV]"):
                line = f.readline()
                i += 1;
                particle.epsilon = eV_to_m*1e-3*mpmathify(line)
            elif (line[:-1] == "# particle sigma [double] [nm]"):
                line = f.readline()
                i += 1;
                particle.sigma = 1e-9*mpmathify(line)
            elif (line[:-1] == "# particle electric charge [double] [e]"):
                line = f.readline()
                i += 1
                particle.charge = e_to_1*mpmathify(line)
            elif (line[:-1] == "# particle start time [double] [s]"):
                line = f.readline()
                i += 1
                startTime1 = s_to_m*mpmathify(line)
            elif (line[:-1] == "# particle start radius [double] [Gm]"):
                line = f.readline()
                i += 1
                startRadius1 = 1e+9*mpmathify(line)
            elif (line[:-1] == "# particle start phi [double] [rad]"):
                line = f.readline()
                i += 1
                startPhi1 = mpmathify(line)
            elif (line[:-1] == "# particle start theta [double] [rad]"):
                line = f.readline()
                i += 1
                startTheta1 = mpmathify(line)
                
            elif (line[:-1] == "# particle r_dot0 [double] [ms^-1 / c]"):
                line = f.readline()
                i += 1
                r_dot0 = mpmathify(line)
            elif (line[:-1] == "# particle r*sintheta*phi_dot0 [double] [ms^-1 / c]"):
                line = f.readline()
                i += 1
                r_sintheta_phi_dot0 = mpmathify(line)
            elif (line[:-1] == "# particle r*theta_dot0 [double] [ms^-1 / c]"):
                line = f.readline()
                i += 1
                r_theta_dot0 = mpmathify(line)
            elif (line[:-1] == "# particle dphi [double] [rad]"):
                line = f.readline()
                i += 1
                dphi = mpmathify(line)
            elif (line[:-1] == "# particle dr/moleculeLength [double]"):
                line = f.readline()
                i += 1
                drCoeff = mpmathify(line)

            # Black hole properties
            elif (line[:-1] == "# Black hole mass [double] [sol]"):
                line = f.readline()
                i += 1
                BH.mass = sol_to_m*mpmathify(line)
            elif (line[:-1] == "# Black hole Kerr parameter (a) [double]"):
                line = f.readline()
                i += 1
                BH.a = mpmathify(line)
            elif (line[:-1] == "# Black hole gravitomagnetic monopole moment (l) [double] [A m]"):
                line = f.readline()
                i += 1
                BH.l = mpmathify(line)
            elif (line[:-1] == "# Black hole electric charge (Q) [double] [e]"):
                line = f.readline()
                i += 1
                BH.charge = e_to_1*mpmathify(line)

            elif (line == "\n" or line==''):
                continue
            else:
                raise Exception("Unrecognised line in properties file:",line)
            
            
            i += 1

            if i >= maxFileSize:
                raise Exception(f'Unsafe file size: number of lines exceeds {maxFileSize}.')
                
    moleculeLength = mpmathify(1.12246204831) * particle.sigma # 2^1/6
                
                
    return BH, particle, moleculeLength