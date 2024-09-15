


    
def writePrimary_cpp(filePath, configureIdentifier:str, motionEquations:list, dofs:tuple,
                    variables:list, enclosers:list, variableLengths:list):
    """
    Write the primary contents (EoMs) to the .cpp file.
    
    Inputs:
        - filePath: Path of the .cpp file.
        - configureIdentifier:str: Indicator of whether the file has been configured.
        - motionEquations:list: EoMs for each dof.
        - dofs:tuple: Names of the degrees of freedom.
        - variables:list: Temporary variables for each dof.
        - enclosers:list: Characters pre and proceeding the temporary variables.
        - variableLengths:list: Number of variables for each dof.
    """
    
    # Write to .cpp file.
    with open(filePath, 'w', encoding='UTF-8') as f:
        f.write(f'{configureIdentifier}')
        f.write("#include \"eulerLagrange.h\"\n\n")

        f.write('const cpp_dec_float_n epsilon_0{ 1. / (4.*M_PI) }; // Geometrised-Gaussian units\n\n\n')
        
        writeIntermediateVariables_cpp(f, variableLengths)
        
        for i in range(len(motionEquations)):
            f.write(f'cpp_dec_float_n calc_{dofs[i]}_next(const BlackHole* BH,\n')
            f.write('\tParticle* p, Particle* ps,\n')
            f.write('\tconst cpp_dec_float_n sigma_p, const cpp_dec_float_n delta_p,\n')
            f.write('\tconst cpp_dec_float_n chi_p, const Precalculated* precalculated,\n')
            f.write('\tconst cpp_dec_float_n dlambda, cpp_dec_float_n* temp)\n')
            f.write('{\n')
            f.write('\t/*\n')
            f.write(f'\tCalculate the next value of {dofs[i]}.\n\n')
            f.write('\tInputs:\n')
            f.write('\t\t- const BlackHole* BH: Instance of the black hole.\n')
            f.write('\t\t- Particle* p: Instance of the particle for which the next coord is being calculated.\n')
            f.write('\t\t- Particle* ps: Instance of the source particle.\n')
            f.write('\t\t- const cpp_dec_float_n sigma_p: Precalculated sigma value.\n')
            f.write('\t\t- const cpp_dec_float_n delta_p: Precalculated delta value.\n')
            f.write('\t\t- const cpp_dec_float_n chi_p: Precalculated chi value.\n')
            f.write('\t\t- const Precalculated* precalculated: Instance of the precalculated variables struct.\n')
            f.write('\t\t- const cpp_dec_float_n dlambda: Lambda step.\n\n')
            
            f.write('\tOutputs:\n')
            f.write(f'\t\t- cpp_dec_float_n {dofs[i]}_new: New coordinate.\n')
            f.write('\t*/\n\n')
            
            
            f.write(f'\tcpp_dec_float_n {dofs[i]}_new' + '{};\n')

            if len(variables[i]) != 0:
                #f.write(f'\tcpp_dec_float_n temp[{len(variables[i])}];\n')
                for j in range(len(variables[i])):
                    f.write(f'\ttemp[{j}] = ' + variables[i][j] + f'; // {enclosers[i][j]}\n')
            f.write('\n\n')

            f.write(f'\t{dofs[i]}_new = 2. * p->{dofs[i]} - p->{dofs[i]}_prev + dlambda*dlambda*({motionEquations[i]});\n')
            f.write(f'\treturn {dofs[i]}_new;\n')
            f.write('}\n\n')

        
        f.write('std::tuple<cpp_dec_float_n, cpp_dec_float_n,\n')
        f.write('\tcpp_dec_float_n, cpp_dec_float_n> eulerMoveMathematica(const BlackHole* BH,\n')
        f.write('\tParticle* p, Particle* ps,\n')
        f.write('\tconst cpp_dec_float_n sigma_p, const cpp_dec_float_n delta_p,\n')
        f.write('\tconst cpp_dec_float_n chi_p, const cpp_dec_float_n dlambda)\n')
        f.write('{\n')
        f.write('\t/*\n')
        f.write(f'\tCalculate the next value of {dofs[i]}.\n\n')
        f.write('\tInputs:\n')
        f.write('\t\t- const BlackHole* BH: Instance of the black hole.\n')
        f.write('\t\t- Particle* p: Instance of the particle for which the next coord is being calculated.\n')
        f.write('\t\t- Particle* ps: Instance of the source particle.\n')
        f.write('\t\t- const cpp_dec_float_n sigma_p: Precalculated sigma value.\n')
        f.write('\t\t- const cpp_dec_float_n delta_p: Precalculated delta value.\n')
        f.write('\t\t- const cpp_dec_float_n chi_p: Precalculated chi value.\n')
        f.write('\t\t- const cpp_dec_float_n dlambda: Lambda step.\n\n')

        f.write('\tOutputs:\n')
        f.write(f'\t\t- cpp_dec_float_n t_new: New t coordinate.\n')
        f.write(f'\t\t- cpp_dec_float_n r_new: New r coordinate.\n')
        f.write(f'\t\t- cpp_dec_float_n phi_new: New phi coordinate.\n')
        f.write(f'\t\t- cpp_dec_float_n theta_new: New theta coordinate.\n')
        f.write('\t*/\n\n')
        f.write('\tcpp_dec_float_n t_new{};\n')
        f.write('\tcpp_dec_float_n r_new{};\n')
        f.write('\tcpp_dec_float_n phi_new{};\n')
        f.write('\tcpp_dec_float_n theta_new{};\n\n')
        
        f.write('\tPrecalculated precalculated{BH, p, ps, sigma_p, delta_p, chi_p, dlambda};\n\n')
        
        
        for i in range(len(motionEquations)):
            f.write(f'\t{dofs[i]}_new = calc_{dofs[i]}_next(BH,\n')
            f.write('\t\tp, ps,\n')
            f.write('\t\tsigma_p, delta_p,\n')
            f.write(f'\t\tchi_p, &precalculated, dlambda')
            if variableLengths[i]!=0:
                f.write(f', precalculated.temp{i}')
            f.write(');\n\n')
                
        
        
        f.write('\treturn {t_new, r_new, phi_new, theta_new};\n')
        f.write('}\n')
    
    
    
    
def writeIntermediateVariables_cpp(f, variableLengths):
    """
    Write intermeditate variables to the .cpp file f.
    
    Inputs:
        - f: Instance of the file.
        - variableLengths:list: Number of temporary variables for each dof.
    """
    
    f.write('Precalculated::Precalculated(const BlackHole* BH,\n')
    f.write('\tParticle* p, Particle* ps,\n')
    f.write('\tconst cpp_dec_float_n sigma_p, const cpp_dec_float_n delta_p,\n')
    f.write('\tconst cpp_dec_float_n chi_p, const cpp_dec_float_n dlambda):\n')
        
    f.write('\tsigma_p2{ sigma_p*sigma_p },\n')
    f.write('\tsigma_p3{ sigma_p*sigma_p2 },\n')
    f.write('\tsigma_p4{ sigma_p2*sigma_p2 },\n')
    f.write('\tdelta_p2{ delta_p*delta_p },\n')
    f.write('\tdelta_p3{ delta_p*delta_p2 },\n')
    f.write('\tdelta_p4{ delta_p2*delta_p2 },\n')
    f.write('\tchi_p2{ chi_p*chi_p },\n')
    f.write('\tchi_p3{ chi_p*chi_p2 },\n')
    f.write('\tchi_p4{ chi_p2*chi_p2 },\n\n')

    f.write('\ttdot{ (p->t_prev - p->t_prevprev) / (dlambda) },\n')
    f.write('\trdot{ (p->r_prev - p->r_prevprev) / (dlambda) },\n')
    f.write('\tphidot{ (p->phi_prev - p->phi_prevprev) / (dlambda) },\n')
    f.write('\tthetadot{ (p->theta_prev - p->theta_prevprev) / (dlambda) },\n\n')
    f.write('\ttsdot{ (ps->t_prev - ps->t_prevprev) / (dlambda) },\n')
    f.write('\trsdot{ (ps->r_prev - ps->r_prevprev) / (dlambda) },\n')
    f.write('\tphisdot{ (ps->phi_prev - ps->phi_prevprev) / (dlambda) },\n')
    f.write('\tthetasdot{ (ps->theta_prev - ps->theta_prevprev) / (dlambda) },\n')
    f.write('\ttsddot{\n')
    f.write('\t\t (ps->t - 2.*ps->t_prev + ps->t_prevprev) / (dlambda*dlambda)},\n')
    f.write('\trsddot{\n')
    f.write('\t\t(ps->r - 2.*ps->r_prev + ps->r_prevprev) / (dlambda*dlambda)},\n')
    f.write('\tphisddot{\n')
    f.write('\t\t(ps->phi - 2.*ps->phi_prev + ps->phi_prevprev) / (dlambda*dlambda)},\n')
    f.write('\tthetasddot{\n')
    f.write('\t\t(ps->theta - 2.*ps->theta_prev + ps->theta_prevprev) / (dlambda*dlambda)},\n\n')

    f.write('\tnx{p->r*sin(p->theta)*cos(p->phi) - ps->r*sin(ps->theta)*cos(ps->phi)},\n') 
    f.write('\tny{p->r*sin(p->theta)*sin(p->phi) - ps->r*sin(ps->theta)*sin(ps->phi)},\n') 
    f.write('\tnz{p->r*cos(p->theta) - ps->r*cos(ps->theta)},\n')
    f.write('\tnmag{sqrt(nx*nx + ny*ny + nz*nz)},\n')
    f.write('\tnmag2{nmag*nmag},\n')
    f.write('\tnmag5{nmag*nmag2*nmag2},\n')
    f.write('\tnmag6{nmag5*nmag},\n')
    f.write('\tnmag7{nmag6*nmag},\n')
    f.write('\tnmag8{nmag6*nmag2},\n')
    f.write('\tnmag12{nmag6*nmag6},\n')
    f.write('\tnmag13{nmag6*nmag7},\n')
    f.write('\tnmag14{nmag7*nmag7},\n\n')


    f.write('\tral_sqr{BH->a2 + BH->l2 + p->r2},\n')
    f.write('\tsintheta{sin(p->theta)},\n')
    f.write('\tsintheta2{sintheta*sintheta},\n')
    f.write('\tsintheta3{sintheta2*sintheta},\n')
    f.write('\tsintheta4{sintheta2*sintheta2},\n')
    f.write('\tsintheta5{sintheta2*sintheta3},\n')
    f.write('\tsintheta6{sintheta3*sintheta3},\n')

    f.write('\t\tsintheta7{sintheta3*sintheta4},\n\n')
    f.write('\tcostheta{cos(p->theta)},\n')
    f.write('\tsin2theta{sin(2.*p->theta)},\n')
    f.write('\tsinphi{sin(p->phi)},\n')
    f.write('\tcosphi{cos(p->phi)},\n\n')

    f.write('\tsinthetas{sin(ps->theta)},\n')
    f.write('\tcosthetas{cos(ps->theta)},\n')
    f.write('\tsinphis{sin(ps->phi)},\n')
    f.write('\tcosphis{cos(ps->phi)},\n')
    f.write('\tsinthetasdot{sin(thetasdot)},\n')
    f.write('\tcosthetasdot{cos(thetasdot)},\n')
    f.write('\tsinphisdot{sin(phisdot)},\n')
    f.write('\tcosphisdot{cos(phisdot)},\n')
    f.write('\tsinthetas2{sinthetas*sinthetas},\n')
    f.write('\tcosthetas2{costhetas*costhetas},\n')
    f.write('\tsinphis2{sinphis*sinphis},\n')
    f.write('\tcosphis2{cosphis*cosphis},\n')
    f.write('\trsdot2{rsdot*rsdot},\n')
    f.write('\tphisdot2{phisdot*phisdot}\n')
    f.write('{\n')
    
    f.write('\t/*\n')
    f.write('\tInitalise the struct.\n')
    f.write('\t*/\n')
    f.write('\treturn;\n')
    f.write('}\n\n\n')
    
    f.write('Precalculated::~Precalculated()\n')
    f.write('{\n')
    
    f.write('\t/*\n')
    f.write('\tDelete the struct.\n')
    f.write('\t*/\n')
    for i in range(4):
        if variableLengths[i] != 0:
            f.write(f'delete[] temp{i};\n')
            f.write(f'temp{i} = nullptr;\n')

    f.write('}\n\n\n')

    
    