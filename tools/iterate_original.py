#!/usr/bin/env python
# Author: Kristjan Haule, March 2007-2017
from scipy import * 
import os,sys,subprocess
"""
This module runs ctqmc impurity solver for one-band model.
The executable should exist in directory params['exe']
"""

Uc=2.25
params = {"exe":   ["./ctqmc",          "# Path to executable"],
          "U":     [Uc,                 "# Coulomb repulsion (F0)"],
          "mu":    [Uc/2.,              "# Chemical potential"],
          "beta":  [100,                "# Inverse temperature"],
          "M" :    [5e6,                "# Number of Monte Carlo steps"],
          "mode":  ["SM",               "# S stands for self-energy sampling, M stands for high frequency moment tail"],
          "cix":   ["one_band.imp",     "# Input file with atomic state"],
          "Delta": ["Delta.inp",        "# Input bath function hybridization"],
          "tsample":[200,               "# how often to record the measurements" ],
          "nom":   [80,                 "# number of Matsubara frequency points to sample"],
        "svd_lmax":[30,                 "# number of SVD functions to project the solution"],
          "aom":   [1,                  "# number of frequency points to determin high frequency tail"],
          "GlobalFlip":[1000000,         "# how often to perform global flip"],
          }

icix="""# Cix file for cluster DMFT with CTQMC
# cluster_size, number of states, number of baths, maximum matrix size
1 4 2 1
# baths, dimension, symmetry, global flip
0       1 0 0
1       1 0 0
# cluster energies for unique baths, eps[k]
0 0
#   N   K   Sz size F^{+,dn}, F^{+,up}, Ea  S
1   0   0    0   1   2         3        0   0
2   1   0 -0.5   1   0         4        0   0.5
3   1   0  0.5   1   4         0        0   0.5
4   2   0    0   1   0         0        0   0
# matrix elements
1  2  1  1    1    # start-state,end-state, dim1, dim2, <2|F^{+,dn}|1>
1  3  1  1    1    # start-state,end-state, dim1, dim2, <3|F^{+,up}|1>
2  0  0  0
2  4  1  1   -1    # start-state,end-state, dim1, dim2, <4|F^{+,up}|2>
3  4  1  1    1
3  0  0  0
4  0  0  0
4  0  0  0
HB2                # Hubbard-I is used to determine high-frequency
# UCoulomb : (m1,s1) (m2,s2) (m3,s2) (m4,s1)  Uc[m1,m2,m3,m4]
0 0 0 0 0.0
# number of operators needed
0
"""

def CreateInputFile(params):
    " Creates input file (PARAMS) for CT-QMC solver"
    f = open('PARAMS', 'w')
    print >> f, '# Input file for continuous time quantum Monte Carlo'
    for p in params:
        print >> f, p, params[p][0], '\t', params[p][1]
    f.close()

def DMFT_SCC(fDelta):
    """This subroutine creates Delta.inp from Gf.out for DMFT on bethe lattice: Delta=t^2*G
    If Gf.out does not exist, it creates Gf.out which corresponds to the non-interacting model
    In the latter case also creates the inpurity cix file, which contains information about
    the atomic states.
    """
    fileGf = 'Gf.out'
    if (os.path.exists(fileGf)): # If output file exists, start from previous iteration
        # Gf = io.read_array(fileGf, columns=(0,-1), lines=(1,-1))
        # In the new Python, io.readarray is dropped and we should use loadtxt instead!
        Gf = loadtxt(fileGf)
    else: # otherwise start from non-interacting limit
        print 'Starting from non-interacting model'
        Gf=[]
        for n in range(2000):
            iom = (2*n+1)*pi/params['beta'][0]
            gf = 2*1j*(iom-sqrt(iom**2+1))
            Gf.append([iom, gf.real, gf.imag])
        Gf = array(Gf)
        # creating impurity cix file
        f = open(params['cix'][0], 'w')
        print >> f, icix
        f.close()
        
    # Preparing input file Delta.inp
    f = open(fDelta, 'w')
    for i in range(len(Gf)):
        print >> f, Gf[i,0], 0.25*Gf[i,1], 0.25*Gf[i,2] # This is DMFT SCC: Delta = t**2*G (with t=1/2)
    f.close()


def Diff(fg1, fg2):
    data1 = loadtxt(fg1).transpose()
    data2 = loadtxt(fg2).transpose()
    diff = sum(abs(data1-data2))/(shape(data1)[0]*shape(data1)[1])
    return diff

# Number of DMFT iterations
Niter = 10

# Creating parameters file PARAMS for qmc execution
CreateInputFile(params)

for it in range(Niter):
    # Constructing bath Delta.inp from Green's function
    DMFT_SCC(params['Delta'][0])

    # Running ctqmc
    print 'Running ---- qmc itt.: ', it, '-----'
    
    subprocess.call(params['exe'][0], shell=True,stdout=sys.stdout,stderr=sys.stderr)
    
    # Some copying to store data obtained so far (at each iteration)
    cmd = 'cp Gf.out Gf.out.'+str(it)
    subprocess.call(cmd, shell=True,stdout=sys.stdout,stderr=sys.stderr)  # copying Gf
    cmd = 'cp Sig.out Sig.out.'+str(it)
    subprocess.call(cmd, shell=True,stdout=sys.stdout,stderr=sys.stderr) # copying Sig
    cmd = 'cp ctqmc.log ctqmc.log.'+str(it)
    subprocess.call(cmd, shell=True,stdout=sys.stdout,stderr=sys.stderr) # copying log file
    
    if it>1:
        diff = Diff('Gf.out', 'Gf.out.'+str(it-1))
        print 'Diff=', diff
        
        
