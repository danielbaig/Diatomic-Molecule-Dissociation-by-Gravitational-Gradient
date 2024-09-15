from numba import njit, types
import re


def simplifyExpressions(motionEquations:'list'):
    """
    Apply the compression algorithm to each of the four dof.
    
    Inputs:
        - motionEquations:list: Expressions for each of the dofs.
        
    Outputs:
        - variables:list: Expressions that have been substituted.
        - motionEquations:list: New expressions with the substitutions.
    """
    
    variables = [[] for i in range(len(motionEquations))]
    enclosers = [[] for i in range(len(motionEquations))]
    
    for i in range(len(motionEquations)):
        print(f'{i+1}/{len(motionEquations)}')
        variables[i], enclosers[i], motionEquations[i] = simplifyExpression(motionEquations[i])
        
    return variables, enclosers, motionEquations


@njit('b1(types.unicode_type)')
def testBrackets(testString:str) -> bool:
    """
    Test to see if all brackets are closed in a string.
    
    Inputs:
        - testString:str: String to test.
        
    Outputs:
        - isClosed:bool: Whether the brackets are closed.
    """
    # Perform checks to see if it is a valid subexpression.
    bracketLevel = 0
    
    for char in testString[1:-1]:
        if char=='(':
            bracketLevel += 1
        elif char==')':
            bracketLevel -= 1
            if bracketLevel == -1:
                return False
    
    return bracketLevel == 0
    
    

def simplifyExpression(expression:str, lowerLength:int=15, upperLength:int=1_000,
                      runNumber:int=1, numTemporary:int=0, tempNameLen:int=9):
    """
    Extract common parts of a C-style expression taking parentheses into account.
    
    Inputs:
        - expression:str: Expression to simplify.
        - lowerLength:int: Lower bound of the length of substring to test.
        - upperLength:int: Upper bound of the length of substring to test.
        - runNumber:int: Number of times the compression procedure has been run.
        - numTemporary:int: Number of temporary variables.
        - tempNameLen:int: Length of the temporary variable.
        
    Outputs:
        - variables:list: Expressions that have been substituted.
        - expression:str: New expression with the substitutions.
    """
    
    
    operators_plus = {'+','-'}
    operators_mul = {'*',}

    stringReplaceThreshold = 200
    numPasses = 30 ##########################
    if upperLength>len(expression)//2: upperLength = len(expression)//2

        
    variables = []
    enclosers = []
    
    if numPasses == 0:
        return variables, expression
    
    
    print(f'\t{runNumber}/{numPasses}')
    i = 0
    
    while i < len(expression) - lowerLength:

        
#         for j in range(0,upperLength - lowerLength):
#             testString = expression[i:i+upperLength-j]
        for j in range(lowerLength, upperLength):
            testString = expression[i:i+j]
                        
            # +...+ / -...- / +...- / -...+
            if (not (testString[0] in operators_plus and testString[-1] in operators_plus 
                and (any(x in testString[1:-1] for x in operators_plus) 
                     or any(x in testString[1:-1] for x in operators_mul))
                and testBrackets(testString)
                #and testString[1]!='>' and expression[i+upperLength-j]!='>')
                and testString[1]!='>' and expression[i+j]!='>')
            # *...*
            and not ((testString[0] in operators_mul or testString[0] in operators_plus) 
                and (testString[-1] in operators_mul or testString[-1] in operators_plus)
                and any(x in testString[1:-1] for x in operators_mul)
                and '+' not in testString[1:-1]
                and testString[1:-1].count('-') <= testString[1:-1].count('>')
                and testBrackets(testString)
                #and testString[1]!='>' and expression[i+upperLength-j]!='>')
                and testString[1]!='>' and expression[i+j]!='>')
            # (...)
           and not (testString[0]=='(' and testString[-1]==')'
                and (any(x in testString[1:-1] for x in operators_plus) 
                     or any(x in testString[1:-1] for x in operators_mul))
                and ',' not in testString
                and expression[i-1] != '/'
                and testBrackets(testString))):
                continue


            
            # Find instances of substring in expression.
            inds = [m.start() for m in re.finditer(re.escape(testString), expression)]
            if len(inds)>1 or (len(testString)<stringReplaceThreshold 
                               and runNumber!=0):
                
                numTemporary += 1

                if numTemporary==11 or numTemporary==101 or numTemporary==1001:
                    tempNameLen += 1
                    assert lowerLength > tempNameLen

                lenTest = len(testString)
                
                
                # Perform substitutions.
                for numSub,k in enumerate(inds):
                    adj = numSub*(tempNameLen - lenTest)
                    
                    if testString[0]=='(' and testString[-1]==')':
                        adj -= 2*numSub
                        expression = (expression[:k + adj] 
                                    + f'temp[{numTemporary-1}]' 
                                      + expression[k + adj + lenTest:])
                    elif testString[0] in operators_plus:
                        expression = (expression[:k + adj] 
                                    + f'+temp[{numTemporary-1}]{testString[-1]}' 
                                      + expression[k + adj + lenTest:])
                    elif testString[0] in operators_mul:
                        expression = (expression[:k + adj] 
                                    + f'*temp[{numTemporary-1}]{testString[-1]}' 
                                      + expression[k + adj + lenTest:])
                
                    
                i += lenTest
                
                

                if testString[0]=='-':
                    variables.append(testString[:-1])
                else:
                    variables.append(testString[1:-1])
                enclosers.append(testString[0] + '...' + testString[-1])

                
                break
                    
        i += 1
                
    if len(variables)==0:
        print(f'No new variables added at run {runNumber}')
            
    elif runNumber < numPasses:
        varNew, enclosersNew, expression = simplifyExpression(expression, runNumber = runNumber + 1,
                                            numTemporary=numTemporary, tempNameLen=tempNameLen)
    
        variables += varNew
        enclosers += enclosersNew

        
    assert len(variables) < 10001, f'Too many temporary variables: {len(variables)}>10001'
        
    return variables, enclosers, expression
                
                
        
        
        
        
def formatExpressions(motionEquations:'list'):
    """
    Replace parts of the expression to make it more efficient and compatable with the already present variables.
    
    Inputs:
        - motionEquations:list: Expressions for each of the dofs.
        
    Outputs:
        - motionEquations:list: New expressions for each of the dofs.
    
    """
    
    patterns = {
        'q' : 'p->charge',
        'sigma(r,theta)': 'sigma_p',
        'delta(r)': 'delta_p',
        'chi(theta)': 'chi_p',
        'Power(Q,2)': 'BH->charge2',
        'Power(l,2)': 'BH->l2',
        'Power(a,2) + Power(l,2) + Power(r,2)': 'precalculated->ral_sqr',
        'Q': 'BH->charge',
        'l': 'BH->l',
        'a': 'BH->a',
        'M': 'BH->mass',
        'Power(a,2)': 'BH->a2',
        'Power(l,2)': 'BH->l2',
        'Power(Q,2)': 'BH->charge2',
        'Power(l,3)': 'BH->l3',
        'Power(a,3)': 'BH->a3',
        'Power(a,4)': 'BH->a4',
        'Power(l,4)': 'BH->l4',
        'Power(a,5)': 'BH->a5',
        'Power(l,5)': 'BH->l5',
        'Power(a,6)': 'BH->a6',
        'Power(l,6)': 'BH->l6',
        'Power(a,7)': 'BH->a7',
        'Power(l,7)': 'BH->l7',
        'Power(a,8)': 'BH->a8',
        'Power(l,8)': 'BH->l8',
        'Power(a,9)': 'BH->a9',
        'Power(l,9)': 'BH->l9',
        'Power(a,10)': 'BH->a10',
        'Cos': 'cos',
        'Sin': 'sin',
        'Csc': '1. / sin',
        'Cot': '1. / tan',
        'Sec': '1. / cos',
        'm' : 'p->mass',
        'r': 'p->r',
        'phi': 'p->phi',
        'theta': 'p->theta',
        'Power(q,2)': '(p->charge*ps->charge)',
        'angular' : 'p->L',
        'energy' : 'p->energy',
        'Power(r,2)': 'p->r2',
        'Power(r,3)': 'p->r3',
        'Power(r,4)': 'p->r4',
        'Power(r,5)': 'p->r5',
        'Power(r,6)': 'p->r6',
        'Power(r,7)': 'p->r7',
        'rs': 'ps->r',
        'phis': 'ps->phi',
        'thetas': 'ps->theta',
        'Power(rdot,2)': '(precalculated->rdot*precalculated->rdot)',
        'Power(phidot,2)': '(precalculated->phidot*precalculated->phidot)',
        'Power(thetadot,2)': '(precalculated->thetadot*precalculated->thetadot)',
        
        'Power(sigma(r,theta),2)' : 'precalculated->sigma_p2',
        'Power(sigma(r,theta),3)' : 'precalculated->sigma_p3',
        'Power(sigma(r,theta),4)' : 'precalculated->sigma_p4',
        'Power(delta(r),2)' : 'precalculated->delta_p2',
        'Power(delta(r),3)' : 'precalculated->delta_p3',
        'Power(delta(r),4)' : 'precalculated->delta_p4',
        'Power(chi(theta),2)' : 'precalculated->chi_p2',
        'Power(chi(theta),3)' : 'precalculated->chi_p3',
        'Power(chi(theta),4)' : 'precalculated->chi_p4',
        
        'Cos(theta)': 'precalculated->costheta',
        'Sin(theta)': 'precalculated->sintheta',
        'Sin(2.*theta)': 'precalculated->sin2theta',
        'Power(Sin(theta),2)': 'precalculated->sintheta2',
        'Power(Sin(theta),3)': 'precalculated->sintheta3',
        'Power(Sin(theta),4)': 'precalculated->sintheta4',
        'Power(Sin(theta),5)': 'precalculated->sintheta5',
        'Power(Sin(theta),6)': 'precalculated->sintheta6',
        'Power(Sin(theta),7)': 'precalculated->sintheta7',
        'Sin(phi)': 'precalculated->sinphi',
        'Cos(phi)': 'precalculated->cosphi',
        
        # Source stuff.
        'Sin(thetas)': 'precalculated->sinthetas',
        'Cos(thetas)': 'precalculated->costhetas',
        'Sin(phis)': 'precalculated->sinphis',
        'Cos(phis)': 'precalculated->cosphis',
        'Sin(phisdot)': 'precalculated->sinphisdot',
        'Cos(phisdot)': 'precalculated->cosphisdot',
        'Sin(thetasdot)': 'precalculated->sinthetasdot',
        'Cos(thetasdot)': 'precalculated->costhetasdot',
        'Power(Sin(thetasdot),2)': 'precalculated->sinthetas2',
        'Power(Cos(thetasdot),2)': 'precalculated->costhetas2',
        'Power(Sin(phisdot),2)': 'precalculated->sinphis2',
        'Power(Cos(phisdot),2)': 'precalculated->cosphis2',
        'Power(rsdot,2)': 'precalculated->rsdot2',
        'Power(phisdot,2)': 'precalculated->phisdot2',
        
        'epsilon': 'ps->epsilon',
        'sigma0': 'ps->sigma',
        'Power(sigma0,6)': 'ps->sigma6',
        'Power(sigma0,12)': 'ps->sigma12',
        
        # n vector
        'Power(nmag,2)': 'precalculated->nmag2',
        'Power(nmag,5)': 'precalculated->nmag5',
        'Power(nmag,6)': 'precalculated->nmag6',
        'Power(nmag,7)': 'precalculated->nmag7',
        'Power(nmag,8)': 'precalculated->nmag8',
        'Power(nmag,12)': 'precalculated->nmag12',
        'Power(nmag,13)': 'precalculated->nmag13',
        'Power(nmag,14)': 'precalculated->nmag14',
        'Power(nx,2)': '(precalculated->nx*precalculated->nx)',
        'Power(ny,2)': '(precalculated->ny*precalculated->ny)',
        'Power(nz,2)': '(precalculated->nz*precalculated->nz)',
        # Keep the same.
        ' 1.*' : '',
        '(1.*' : '(',
        '2.' : '2.',
        '/2' : '/2.',
        'Power' : 'pow',
        'epsilon0' : 'epsilon_0',
        'nmag' : 'precalculated->nmag',
        'tdot': 'precalculated->tdot',
        'rdot' : 'precalculated->rdot',
        'phidot' : 'precalculated->phidot',
        'thetadot' : 'precalculated->thetadot',
        'tsdot': 'precalculated->tsdot',
        'rsdot' : 'precalculated->rsdot',
        'phisdot' : 'precalculated->phisdot',
        'thetasdot' : 'precalculated->thetasdot',
        'tsddot':'precalculated->tsddot',
        'rsddot' : 'precalculated->rsddot',
        'phisddot':'precalculated->phisddot',
        'thetasddot' : 'precalculated->thetasddot',
        'nx' : 'precalculated->nx',
        'ny' : 'precalculated->ny',
        'nz' : 'precalculated->nz',
        
        # Plane case.
        'sigma(r,Pi/2.)': 'sigma_p',
        'chi(Pi/2.)': 'chi_p',
        'Power(sigma(r,Pi/2.),2)' : 'precalculated->sigma_p2',
        'Power(sigma(r,Pi/2.),3)' : 'precalculated->sigma_p3',
        'Power(sigma(r,Pi/2.),4)' : 'precalculated->sigma_p4',
        'Power(chi(Pi/2.),2)' : 'precalculated->chi_p2',
        'Power(chi(Pi/2.),3)' : 'precalculated->chi_p3',
        'Power(chi(Pi/2.),4)' : 'precalculated->chi_p4',
    }
    
    maxPatternLength = max(map(len, patterns.keys()))
    
    for i in range(len(motionEquations)):
        j = 0
        while j < len(motionEquations[i]):
            # Run from largest pattern lengt to smallest.
            for k in range(0,maxPatternLength):
                section = motionEquations[i][j:j+maxPatternLength - k] 
                
                # Replace part.
                if section in patterns.keys():
                    motionEquations[i] = (motionEquations[i][:j] + patterns[section] 
                                        + motionEquations[i][j+maxPatternLength - k:])
                    j += len(patterns[section]) - 1
                    break
                    
            j += 1
            
        # Remove spaces.
        motionEquations[i] = motionEquations[i].replace(" ", "")
                
                
                    
    return motionEquations