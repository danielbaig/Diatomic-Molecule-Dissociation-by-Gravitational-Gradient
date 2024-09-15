import numpy as np
from numpy import linalg

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from tqdm import tqdm

from mpmath import mp, mpf, sin, cos, sqrt, acos

from pathlib import Path
pathtohere = Path.cwd()

epsilon_0 = 1. / (4.*np.pi) # Geometrised-Gaussian units




def displaySystem(BH, coords1:np.ndarray, coords2:np.ndarray,
                  isPrecise:bool=True, renderFull:bool=True):
    """
    Display the particle(')s(') motion through the system.
    
    Inputs:
        - BH: Instance of the black hole.
        - coords1:np.ndarray: History of the coordinates of particle 1.
        - coords2:np.ndarray: History of the coordinates of particle 2.
        - isPrecise:bool: Whether to use the ultra precise measurements.
    """
    
    r_Q2 = BH.charge*BH.charge / (4*np.pi*epsilon_0)
    r_s = 2.*BH.mass
    
    fig = plt.figure(figsize=(6,6), dpi=150)
    ax = fig.add_subplot(111, projection='3d')
    
    # Sphere mesh.
    phi = np.linspace(0, 2 * np.pi, 100)
    theta = np.linspace(0, np.pi, 100)
    x = np.cos(phi)[:, None] * np.sin(theta)[None,:]
    y = np.sin(phi)[:,None] * np.sin(theta)[None,:]
    z = np.ones(np.size(phi))[:,None] * np.cos(theta)[None,:]
    
    # Inner horizon
    r = r_s / 2. - np.sqrt(r_s*r_s / 4. - BH.a*BH.a - r_Q2)
    ax.plot_surface(r * x, r * y, r * z, color='black', alpha=0.2)
    print(f'Inner horizon: {r}')
    # Outer horizon
    r = r_s / 2. + np.sqrt(r_s*r_s / 4. - BH.a*BH.a - r_Q2)
    ax.plot_surface(r * x, r * y, r * z, color='m', alpha=0.2)
    print(f'Outer horizon: {r}')
    
    # Inner ergosphere
    r = r_s / 2. - np.sqrt(r_s*r_s / 4. - BH.a*BH.a*np.cos(theta)*np.cos(theta) - r_Q2)
    ax.plot_surface(r[None,:] * x, r[None,:] * y, r[None,:] * z, color='cyan', alpha=0.2)
    print(f'Inner ergosphere: [{r.min()},{r.max()}]')
    # Outer ergosphere
    r = r_s / 2. + np.sqrt(r_s*r_s / 4. - BH.a*BH.a*np.cos(theta)*np.cos(theta) - r_Q2)
    ax.plot_surface(r[None,:] * x, r[None,:] * y, r[None,:] * z, color='yellow', alpha=0.2)
    print(f'Outer ergosphere: [{r.min()},{r.max()}]')
    
    if not isPrecise:
        # Particle1 trajectory.
        x1 = coords1[:,2] * np.sin(coords1[:,4]) * np.cos(coords1[:,3])
        y1 = coords1[:,2] * np.sin(coords1[:,4]) * np.sin(coords1[:,3])
        z1 = coords1[:,2] * np.cos(coords1[:,4])

        # Particle2 trajectory.
        x2 = coords2[:,2] * np.sin(coords2[:,4]) * np.cos(coords2[:,3])
        y2 = coords2[:,2] * np.sin(coords2[:,4]) * np.sin(coords2[:,3])
        z2 = coords2[:,2] * np.cos(coords2[:,4])
        
    else:
        x1 = np.array([c2 * sin(c4) * cos(c3) 
                      for c2, c4, c3 in zip(coords1[:,2], coords1[:,4], coords1[:,3])], dtype=object)
        y1 = np.array([c2 * sin(c4) * sin(c3) 
                      for c2, c4, c3 in zip(coords1[:,2], coords1[:,4], coords1[:,3])], dtype=object)
        z1 = np.array([c2 * cos(c4) 
                      for c2, c4 in zip(coords1[:,2], coords1[:,4])], dtype=object)
        
        x2 = np.array([c2 * sin(c4) * cos(c3) 
                      for c2, c4, c3 in zip(coords2[:,2], coords2[:,4], coords2[:,3])], dtype=object)
        y2 = np.array([c2 * sin(c4) * sin(c3) 
                      for c2, c4, c3 in zip(coords2[:,2], coords2[:,4], coords2[:,3])], dtype=object)
        z2 = np.array([c2 * cos(c4) 
                      for c2, c4 in zip(coords2[:,2], coords2[:,4])], dtype=object)
        
    ax.scatter(x1,y1,z1, c='b', marker='.',s=1.)
    
    if renderFull:
        ax.scatter(x2,y2,z2, c='g', marker='.',s=1.)
    
    
    # Create appropiate labels.
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    
    # Limit graph size.
    size = float(1.1*coords1[0,2])
    ax.set_xlim(-size, size)
    ax.set_ylim(-size, size)
    ax.set_zlim(-size, size)
    
    #plt.show()
    plt.savefig(pathtohere / 'plots/system3D.png', bbox_inches='tight')
    plt.close(fig)
    
def displayMolecule(coords1:np.ndarray, coords2:np.ndarray, isPrecise:bool=True):
    """
    Display the planar movement of the atoms in the molecule.
    
    Inputs:
        - coords1:np.ndarray: History of the coordinates of particle 1.
        - coords2:np.ndarray: History of the coordinates of particle 2.
        - isPrecise:bool: Whether to use the ultra precise measurements.
    """
    
    fig = plt.figure(figsize=(6,6), dpi=150)
    ax = fig.add_subplot()
    
    if not isPrecise:
        # Particle1 trajectory.
        x1 = coords1[:,2] * np.sin(coords1[:,4]) * np.cos(coords1[:,3])
        y1 = coords1[:,2] * np.sin(coords1[:,4]) * np.sin(coords1[:,3])

        # Particle2 trajectory.
        x2 = coords2[:,2] * np.sin(coords2[:,4]) * np.cos(coords2[:,3])
        y2 = coords2[:,2] * np.sin(coords2[:,4]) * np.sin(coords2[:,3])
        
    else:
        x1 = np.array([c2 * sin(c4) * cos(c3) 
                      for c2, c4, c3 in zip(coords1[:,2], coords1[:,4], coords1[:,3])], dtype=object)
        y1 = np.array([c2 * sin(c4) * sin(c3) 
                      for c2, c4, c3 in zip(coords1[:,2], coords1[:,4], coords1[:,3])], dtype=object)
        
        x2 = np.array([c2 * sin(c4) * cos(c3) 
                      for c2, c4, c3 in zip(coords2[:,2], coords2[:,4], coords2[:,3])], dtype=object)
        y2 = np.array([c2 * sin(c4) * sin(c3) 
                      for c2, c4, c3 in zip(coords2[:,2], coords2[:,4], coords2[:,3])], dtype=object)
   
    COM = np.asarray([x1+x2, y1+y2]) / 2.

    xoffset = min(min(x1),min(x2))

    ax.scatter(x1 - COM[0],y1 - COM[1], c='b', marker='.',s=1., zorder=3)
    ax.scatter(x2 - COM[0],y2 - COM[1], c='g', marker='.',s=1., zorder=3)
    
    # Draw bond.
    for i in tqdm(range(len(x1))):
        if i%(len(x1)//1000) == 0 and i>0.67*len(x1) and i < 0.69*len(x1):
            ax.plot((x1[i] - COM[0,i], x2[i] - COM[0,i]),(y1[i] - COM[1,i],y2[i] - COM[1,i]),
                    c='k', alpha=0.9*((i/len(x1)-0.67)/(0.69-0.67)), zorder=2)
        elif i==0:
            ax.plot((x1[i] - COM[0,i], x2[i] - COM[0,i]),(y1[i] - COM[1,i], y2[i] - COM[1,i]), c='r', alpha=0.9, zorder=4)

    
    
    # Create appropiate labels.
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title('CoM frame')
    
    #plt.show()
    plt.savefig(pathtohere / 'plots/molecule.png', bbox_inches='tight')
    plt.close(fig)
    
    
def displayCoordinateStats(BH, coords1:np.ndarray, coords2:np.ndarray, isPrecise:bool=True):
    """
    Display how each the coordinates changes for both particles.
   
    Inputs:
        - BH: Instance of the black hole.
        - coords1:np.ndarray: History of the coordinates of particle 1.
        - coords2:np.ndarray: History of the coordinates of particle 2.
        - isPrecise:bool: Whether to use the ultra precise measurements.
    """
    
    r_Q2 = BH.charge*BH.charge / (4*np.pi*epsilon_0)
    r_s = 2.*BH.mass
    
    lambdas = coords1[:,0]
    coords = (coords1[:,1:] + coords2[:,1:]) / 2.
    
    colours = ('r','g','b','m')
    ylabels = ('t','r',r'$\phi$',r'$\theta$')
    
    
    fig = plt.figure(figsize=(8,8), tight_layout=True)
    
    for i in range(4):
        ax = fig.add_subplot(2,2,i+1)
        offset = min(coords[:,i]) if i==1 else 0
        
        ax.scatter(lambdas, coords[:,i] - offset, marker='.',c=colours[i])
        if i==1:
            if not isPrecise:
                r_erg = r_s / 2. + np.sqrt(r_s*r_s / 4. - BH.a*BH.a*np.cos(coords[:,3])*np.cos(coords[:,3]) - r_Q2)
            else:
                r_erg = r_s / 2. + np.array([sqrt(r_s*r_s / 4. - BH.a*BH.a*cos(c3)*cos(c3) - r_Q2)
                      for c3 in coords[:,3]], dtype=object)
            
            r_eventHorizon = r_s / 2. + np.sqrt(r_s*r_s / 4. - BH.a*BH.a - r_Q2)
            
            if min(coords[:,i]) <= max(r_erg):
                ax.plot(lambdas,r_erg-offset,c='y')
                ax.axhline(r_eventHorizon-offset,c='m')
            else:
                ax.set_ylim(float(min(coords[:,i])-offset), float(max(coords[:,i])-offset))
                ax.set_title(f'+{offset}')
        
        ax.grid()
        
        # Create appropiate labels.
        ax.set_xlabel(r'$\lambda$')
        ax.set_ylabel(ylabels[i])
        
    #plt.show()
    plt.savefig(pathtohere / 'plots/coordStats.png', bbox_inches='tight')
    plt.close(fig)
    
    
def displayPhaseSpace_phi_r(BH, coords1:np.ndarray, coords2:np.ndarray, isPrecise:bool=True):
    """
    Display how the particles move in the phi-r phase space.
   
    Inputs:
        - BH: Instance of the black hole.
        - coords1:np.ndarray: History of the coordinates of particle 1.
        - coords2:np.ndarray: History of the coordinates of particle 2.
        - isPrecise:bool: Whether to use the ultra precise measurements.
    """
    
    r_Q2 = BH.charge*BH.charge / (4*np.pi*epsilon_0)
    r_s = 2.*BH.mass    

    
    # Create figure.
    fig = plt.figure(figsize=(8,8))
    
    ax = fig.add_subplot()
    offset = min(min(coords1[:,2]),min(coords2[:,2]))

    ax.scatter(coords1[:,3], coords1[:,2] - offset, marker='.',c='b')
    ax.scatter(coords2[:,3], coords2[:,2] - offset, marker='.',c='g')
    
    if not isPrecise:
        r_erg = r_s / 2. + np.sqrt(r_s*r_s / 4. - BH.a*BH.a*np.cos(coords1[:,4])*np.cos(coords1[:,4]) - r_Q2)
    else:
        r_erg = r_s / 2. + np.array([sqrt(r_s*r_s / 4. - BH.a*BH.a*cos(c3)*cos(c3) - r_Q2)
              for c3 in coords1[:,4]], dtype=object)

    r_eventHorizon = r_s / 2. + np.sqrt(r_s*r_s / 4. - BH.a*BH.a - r_Q2)

    if offset <= max(r_erg):
        ax.plot(coords1[:,3],r_erg-offset,c='y')
        ax.axhline(r_eventHorizon-offset,c='m')
    else:
        ax.set_ylim(0, float(max(max(coords1[:,2]),max(coords2[:,2])) - offset))
        ax.set_title(f'+{offset}')
        
    ax.grid()

    # Create appropiate labels.
    ax.set_xlabel(r'$\phi$')
    ax.set_ylabel('r')
        
    #plt.show()
    plt.savefig(pathtohere / 'plots/phaseSpace_phi_r.png', bbox_inches='tight')
    plt.close(fig)
    
    
def displayMoleculeStats(coords1:np.ndarray, coords2:np.ndarray,
                         moleculeLength:float, isPrecise:bool=True):
    """
    Display how the moleculear angle and length vary.
    
    
    Inputs:
        - coords1:np.ndarray: History of the coordinates of particle 1.
        - coords2:np.ndarray: History of the coordinates of particle 2.
        - moleculeLength:float: The average molecule length.
        - isPrecise:bool: Whether to use the ultra precise measurements.
    """
    
    
    if not isPrecise:
        # Particle1 trajectory.
        x1 = coords1[:,2] * np.sin(coords1[:,4]) * np.cos(coords1[:,3])
        y1 = coords1[:,2] * np.sin(coords1[:,4]) * np.sin(coords1[:,3])
        z1 = coords1[:,2] * np.cos(coords1[:,4])

        # Particle2 trajectory.
        x2 = coords2[:,2] * np.sin(coords2[:,4]) * np.cos(coords2[:,3])
        y2 = coords2[:,2] * np.sin(coords2[:,4]) * np.sin(coords2[:,3])
        z2 = coords2[:,2] * np.cos(coords2[:,4])
        
    else:
        x1 = np.array([c2 * sin(c4) * cos(c3) 
                      for c2, c4, c3 in zip(coords1[:,2], coords1[:,4], coords1[:,3])], dtype=object)
        y1 = np.array([c2 * sin(c4) * sin(c3) 
                      for c2, c4, c3 in zip(coords1[:,2], coords1[:,4], coords1[:,3])], dtype=object)
        z1 = np.array([c2 * cos(c4) 
                      for c2, c4 in zip(coords1[:,2], coords1[:,4])], dtype=object)
        
        x2 = np.array([c2 * sin(c4) * cos(c3) 
                      for c2, c4, c3 in zip(coords2[:,2], coords2[:,4], coords2[:,3])], dtype=object)
        y2 = np.array([c2 * sin(c4) * sin(c3) 
                      for c2, c4, c3 in zip(coords2[:,2], coords2[:,4], coords2[:,3])], dtype=object)
        z2 = np.array([c2 * cos(c4) 
                      for c2, c4 in zip(coords2[:,2], coords2[:,4])], dtype=object)
    
    COM = np.asarray([x1+x2, y1+y2, z1+z2]) / 2.
    delta = np.asarray([x2-x1, y2-y1, z2-z1])
    
    
    if not isPrecise:
        angle = np.arccos((COM[0]*delta[0] + COM[1]*delta[1] + COM[2]*delta[2]) 
                          / linalg.norm(COM, axis=0) / linalg.norm(delta, axis=0))
        
    else:
        # Compute the dot product of COM and delta
        dot_product = sum(c * d for c, d in zip(COM, delta))
        
        
        # Compute the norms of COM and delta
        norm_COM = np.array([sqrt(s) for s in sum(c*c for c in COM)], dtype=object)
        norm_delta = np.array([sqrt(s) for s in sum(d*d for d in delta)], dtype=object)
        
        
        # Compute the cosine argument
        cos_argument = dot_product / (norm_COM * norm_delta)

        # Compute the arccosine using mpmath.acos
        # Ensure the value is within the valid range for acos
        cos_argument = np.array([max(-1, m) for m in [min(1, n) for n in cos_argument]],
                                dtype=object)  # Clamp value to [-1, 1]

        # Calculate the angle in radians
        angle = np.array([acos(a) for a in cos_argument], dtype=object)
    
    
    
    # Generate figure.
    fig = plt.figure(figsize=(8,10), tight_layout=True)
    
    
    # Angular.
    ax = fig.add_subplot(2,1,1)
    ax.scatter(coords1[:,0], angle / (2.*np.pi), c='c', marker='.')
    ax.grid()
    
    #ax.set_ylim(0, np.pi)
    ax.set_xlabel(r'$\lambda$')
    ax.set_ylabel(r'angle out of radial axis [$\tau$rad]')
    
    # Separation.
    ax = fig.add_subplot(2,1,2)
    if not isPrecise:
        ax.scatter(coords1[:,0], linalg.norm(delta,axis=0), c='lime', marker='.')
    else:
        separationOffset = min(norm_delta)
        ax.scatter(coords1[:,0], norm_delta - separationOffset, c='lime', marker='.')
    ax.grid()
    
    # Only include equilibrium line if it is passed.
    if any(norm_delta < moleculeLength) and any(norm_delta>moleculeLength):
        ax.axhline(moleculeLength, c='r')
    
    ax.set_xlabel(r'$\lambda$')
    ax.set_ylabel(f'molecule length\n+{separationOffset}')
    
    #plt.show()
    plt.savefig(pathtohere / 'plots/moleculeStats.png', bbox_inches='tight')
    plt.close(fig)
        