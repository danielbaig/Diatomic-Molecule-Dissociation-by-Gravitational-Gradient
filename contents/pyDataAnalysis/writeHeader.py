
decimalPlaces = 40

def writePrimary_h(filePath, motionEquations:list, dofs:tuple,
                    variables:list, variableLengths:list):
    """
    Write the primary contents (EoMs) to the .cpp file.
    
    Inputs:
        - filePath: Path of the .h file.
        - motionEquations:list: EoMs for each dof.
        - dofs:tuple: Names of the degrees of freedom.
        - variables:list: Temporary variables for each dof.
        - variableLengths:list: Number of variables for each dof.
    """
    
    with open(filePath, 'w', encoding='UTF-8') as f:
        f.write('#ifndef EULERLAGRANGE_H\n')
        f.write('#define EULERLAGRANGE_H\n\n')
        
        f.write('#include <iostream>\n')
        f.write('#include <tuple>\n')
        f.write('#define _USE_MATH_DEFINES\n')
        f.write('#include "math.h"\n\n')
        f.write('#include <boost/multiprecision/cpp_dec_float.hpp>\n')
        f.write('using namespace boost::multiprecision;\n')
        f.write(f'using cpp_dec_float_n = number<cpp_dec_float<{decimalPlaces}>>;\n\n')


        f.write('#include "blackHole.h"\n')
        f.write('#include "particle.h"\n\n\n')
        
        
        writeIntermediateVariables_h(f, variableLengths)
        
        for i in range(len(motionEquations)):
            f.write(f'cpp_dec_float_n calc_{dofs[i]}_next(const BlackHole* BH,\n')
            f.write('\tParticle* p, Particle* ps,\n')
            f.write('\tconst cpp_dec_float_n sigma_p, const cpp_dec_float_n delta_p,\n')
            f.write('\tconst cpp_dec_float_n chi_p, const Precalculated* precalculated,\n')
            f.write('\tconst cpp_dec_float_n dlambda,\n')
            f.write('\tcpp_dec_float_n* temp=nullptr);\n\n')
        
        f.write('std::tuple<cpp_dec_float_n, cpp_dec_float_n,')
        f.write('\tcpp_dec_float_n, cpp_dec_float_n> eulerMoveMathematica(const BlackHole* BH,\n')
        f.write('\tParticle* p, Particle* ps,\n')
        f.write('\tconst cpp_dec_float_n sigma_p, const cpp_dec_float_n delta_p,')
        f.write('\tconst cpp_dec_float_n chi_p, const cpp_dec_float_n dlambda);\n\n')

        f.write('#endif EULERLAGRANGE_H\n')
    


def writeIntermediateVariables_h(f, variableLengths):
    """
    Write intermeditate variables to the header (.h) file f.
    
    Inputs:
        - f: Instance of the file.
        - variableLengths:list: Number of temporary variables for each dof.
    """

    f.write('struct Precalculated{\n')
    
    for i in range(4):
        if variableLengths[i] != 0:
            f.write(f'\tcpp_dec_float_n* temp{i}' + '{ new cpp_dec_float_n[' + f'{variableLengths[i]}' + ']{} };\n')
    f.write('\n\n')
    
    f.write('\tconst cpp_dec_float_n sigma_p2{};\n')
    f.write('\tconst cpp_dec_float_n sigma_p3{};\n')
    f.write('\tconst cpp_dec_float_n sigma_p4{};\n')
    f.write('\tconst cpp_dec_float_n delta_p2{};\n')
    f.write('\tconst cpp_dec_float_n delta_p3{};\n')
    f.write('\tconst cpp_dec_float_n delta_p4{};\n')
    f.write('\tconst cpp_dec_float_n chi_p2{};\n')
    f.write('\tconst cpp_dec_float_n chi_p3{};\n')
    f.write('\tconst cpp_dec_float_n chi_p4{};\n\n')

    f.write('\tconst cpp_dec_float_n tdot{};\n')
    f.write('\tconst cpp_dec_float_n rdot{};\n')
    f.write('\tconst cpp_dec_float_n phidot{};\n')
    f.write('\tconst cpp_dec_float_n thetadot{};\n\n')
    f.write('\tconst cpp_dec_float_n tsdot{};\n')
    f.write('\tconst cpp_dec_float_n rsdot{};\n')
    f.write('\tconst cpp_dec_float_n phisdot{};\n')
    f.write('\tconst cpp_dec_float_n thetasdot{};\n')
    f.write('\tconst cpp_dec_float_n tsddot{};\n')
    f.write('\tconst cpp_dec_float_n rsddot{};\n')
    f.write('\tconst cpp_dec_float_n phisddot{};\n')
    f.write('\tconst cpp_dec_float_n thetasddot{};\n\n')

    f.write('\tconst cpp_dec_float_n nx{};\n') 
    f.write('\tconst cpp_dec_float_n ny{};\n') 
    f.write('\tconst cpp_dec_float_n nz{};\n')
    f.write('\tconst cpp_dec_float_n nmag{};\n')
    f.write('\tconst cpp_dec_float_n nmag2{};\n')
    f.write('\tconst cpp_dec_float_n nmag5{};\n')
    f.write('\tconst cpp_dec_float_n nmag6{};\n')
    f.write('\tconst cpp_dec_float_n nmag7{};\n')
    f.write('\tconst cpp_dec_float_n nmag8{};\n')
    f.write('\tconst cpp_dec_float_n nmag12{};\n')
    f.write('\tconst cpp_dec_float_n nmag13{};\n')
    f.write('\tconst cpp_dec_float_n nmag14{};\n\n')


    f.write('\tconst cpp_dec_float_n ral_sqr{};\n')
    f.write('\tconst cpp_dec_float_n sintheta{};\n')
    f.write('\tconst cpp_dec_float_n sintheta2{};\n')
    f.write('\tconst cpp_dec_float_n sintheta3{};\n')
    f.write('\tconst cpp_dec_float_n sintheta4{};\n')
    f.write('\tconst cpp_dec_float_n sintheta5{};\n')
    f.write('\tconst cpp_dec_float_n sintheta6{};\n')
    f.write('\tconst cpp_dec_float_n sintheta7{};\n\n')

    f.write('\tconst cpp_dec_float_n costheta{};\n')
    f.write('\tconst cpp_dec_float_n sin2theta{};\n')
    f.write('\tconst cpp_dec_float_n sinphi{};\n')
    f.write('\tconst cpp_dec_float_n cosphi{};\n\n')

    f.write('\tconst cpp_dec_float_n sinthetas{};\n')
    f.write('\tconst cpp_dec_float_n costhetas{};\n')
    f.write('\tconst cpp_dec_float_n sinphis{};\n')
    f.write('\tconst cpp_dec_float_n cosphis{};\n')
    f.write('\tconst cpp_dec_float_n sinthetasdot{};\n')
    f.write('\tconst cpp_dec_float_n costhetasdot{};\n')
    f.write('\tconst cpp_dec_float_n sinphisdot{};\n')
    f.write('\tconst cpp_dec_float_n cosphisdot{};\n')
    f.write('\tconst cpp_dec_float_n sinthetas2{};\n')
    f.write('\tconst cpp_dec_float_n costhetas2{};\n')
    f.write('\tconst cpp_dec_float_n sinphis2{};\n')
    f.write('\tconst cpp_dec_float_n cosphis2{};\n')
    f.write('\tconst cpp_dec_float_n rsdot2{};\n')
    f.write('\tconst cpp_dec_float_n phisdot2{};\n\n')
    
    f.write('Precalculated(const BlackHole* BH,\n')
    f.write('\t\tParticle* p, Particle* ps,\n')
    f.write('\t\tconst cpp_dec_float_n sigma_p, const cpp_dec_float_n delta_p,\n')
    f.write('\t\tconst cpp_dec_float_n chi_p, const cpp_dec_float_n dlambda);\n\n')
    
    f.write('~Precalculated();\n')
    
    f.write('};\n\n')
    