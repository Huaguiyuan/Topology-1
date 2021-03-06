# this is working fine. we have produced result for t20.0 using this. it give the same transition U as triqs code.
#!/usr/bin/env python
# Author: soumen, 13/03/16
from numpy import * 
import os,sys,subprocess
import numpy as np
import utility_IHM_dos_haule
Uc = 0.0
beta = 0.0
M = 10e6
#/home/soumen/INSTALL-PACKAGES/dmft_w2k_2015/ctqmc'
#/home/soumen/INSTALLPACKAGES/mpich2/bin/mpirun
params = {"exe"          : ['./ctqmc'          , "# Path to executable"],
          "Delta"        : ["Delta.inp"         , "# Input bath function hybridization"],
          "cix"          : ["one_band.imp"      , "# Input file with atomic state"],
          "U"            : [Uc                  , "# Coulomb repulsion (F0)"],
          "mu"           : [Uc/2               , "# Chemical potential"],
          "beta"         : [beta                , "# Inverse temperature"],
          "M"            : [M                  , "# Number of Monte Carlo steps"],
          "nom"          : [180                , "# number of Matsubara frequency points to sample"],
          "aom"          : [1                   , "# number of frequency points to determin high frequency tail"],
          "tsample"      : [200                  , "# how often to record the measurements" ],
          "mode"	 : ["SM"		,"# S stands for self-energy sampling, M stands for high frequency moment tail"],
          "svd_lmax"     : [30                  ,"# number of SVD functions to project the solution"]
}



icix="""# Cix file for cluster DMFT with CTQMC
# cluster_size, number of states, number of baths, maximum matrix size
1 4 2 1
# baths, dimension, symmetry, global flip
0       1 0 0
1       1 1 0
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

def mix_deltaA(ita,mix):
	dold = loadtxt('DeltaA.inp.'+str(ita-1))
	dnew = loadtxt('Delta.inp')
	delta = open('Delta.inp','w')
	data = mix*dnew + (1-mix)*dold
	for i in range(N):
		print >> delta, dold[i][0], data[i][1], data[i][2], data[i][3], data[i][4]
	delta.close()
	
def mix_deltaB(ita,mix):
	dold = loadtxt('DeltaB.inp.'+str(ita-1))
	dnew = loadtxt('Delta.inp')
	delta = open('Delta.inp','w')
	data = mix*dnew + (1-mix)*dold
	for i in range(N):
		print >> delta, dold[i][0], data[i][1], data[i][2], data[i][3], data[i][4]
	delta.close()


def CreateInputFileA(params):
	" Creates input file (PARAMS) for CT-QMC solver"
	fA = open('PARAMSA', 'w')
	f = open('PARAMS', 'w')
	print >> fA, '# Input file for continuous time quantum Monte Carlo'
	print >> f, '# Input file for continuous time quantum Monte Carlo'
	for p in params:
		print >> f, p, params[p][0], '\t', params[p][1]
		print >> fA, p, params[p][0], '\t', params[p][1]
	f.close()

def CreateInputFileB(params):
	    " Creates input file (PARAMS) for CT-QMC solver"
	    f = open('PARAMS', 'w')
	    fB = open('PARAMSB', 'w')
	    print >> f, '# Input file for continuous time quantum Monte Carlo'
	    print >> fB, '# Input file for continuous time quantum Monte Carlo'
	    for p in params:
		print >> f, p, params[p][0], '\t', params[p][1]
		print >> fB, p, params[p][0], '\t', params[p][1]
	    f.close()

def create_cix_input(params,icix):
		f = open(params['cix'][0], 'w')
		print >> f, icix
		f.close()
		

#setting the environment veriables
mpirun_desktop = '/home/soumen/INSTALLPACKAGES/mpich2/bin/mpirun'
ctqmc_desktop  = '/home/soumen/INSTALLPACKAGES/dmft_w2k/ctqmc'
mpirun_rahman  = '/home/soumen/INSTALL-PACKAGES/mpich2/bin/mpirun'
ctqmc_rahman   = '/home/soumen/INSTALL-PACKAGES/dmft_w2k_2015/ctqmc'
#params['exe'][0] = ctqmc_rahman
mpi = mpirun_rahman
# following mpi for ./hybri
mpi2 = '/home/soumen/softwares/mpich-3.1.4/bin/mpirun' 


#intialising the green function
params['M'][0] = 10e6
beta=100.0
delta = 0.5
dos_file = "square.dat"
Niter = 20
N = 1025                     # no of matsubara frequency
params['beta'][0]=beta
g = utility_IHM_dos_haule.green(beta,delta)
#g.non_iteracting_G_IHM(dos_file, 0.0)
'''
it_previous = 59 
cmd = 'cp SigA.out.'+str(it_previous) +  ' SigA.out.0'
print os.popen(cmd).read()
cmd = 'cp SigB.out.'+str(it_previous) +  ' SigB.out.0'
print os.popen(cmd).read()
'''
g.create_cix_input(params,icix)
U_list=[0.5]
#it_previous = 59
for i,U in enumerate(U_list):

	params['U'][0]=U 
        fh_info = open('info.dat', 'w')
        hyb_info = open('hyb_info.dat','w')
	directory = 'U%s'%U
        cmd = 'mkdir %s'%directory
        print os.popen(cmd).read()
	chemical_potential = U/2.
	#it_list = range(it_previous+1,Niter)
        #if U==0.0: Niter = 10
        #else: Niter = 50
        #if U>0.5: Niter = 60
	#if i==0: it_list = range(it_previous+1,Niter)
	#else: it_list = range(Niter)
	#l = 0

	
        #if previous U self energy is there
        if i!=0 :
		it_previous = Niter-1
		cmd = 'cp SigA.out.'+str(it_previous) +  ' SigA.out.0'
		print os.popen(cmd).read()
		cmd = 'cp SigB.out.'+str(it_previous) +  ' SigB.out.0'
		print os.popen(cmd).read()

        #if U>0.5: Niter = 60
        #if i==0: it_list = range(it_previous+1,Niter)
        it_list = range(Niter)
        l = 0
	for it in it_list:

		#chemical potential update details
		update_start = 60 
		update_step = 7
		update_end = 9 #update will stop after Niter-update_end
		if it<70: coeffi = 0.3
                elif it< 120: coeffi = 0.2
                else: coeffi = 0.1
                if it<100 : mix = 0.7
                else : mix = 0.4
                l = l+1
		
		# Constructing Delta.inp for A sublattice 
		params['mu'][0]=chemical_potential - delta
		CreateInputFileA(params)
                cmd = "aprun -j 1 -n 24 -N 24 " +  "./hybri A " + str(chemical_potential) + " "+str(it)
                subprocess.call(cmd,shell=True,stdout=hyb_info,stderr=hyb_info)
   		hyb_info.flush()
    		if it !=0: mix_deltaA(it,mix)
		
		# Running ctqmc for A
    		print 'Running ---- qmc itt.A: ', it, '--U,coeffi,mix---',U,coeffi,mix
                #aprun -j 1 -n 24 -N 24 ./ctqmc PARAMS >& aprun_utput.txt
   		cmd = "aprun -j 1 -n 24 -N 24 "+params['exe'][0]+'  PARAMS > nohup_imp.out 2>&1 '
    		subprocess.call(cmd,shell=True,stdout=fh_info,stderr=fh_info)
   		fh_info.flush()

    		#saving A
		cmd = 'cp Delta.inp DeltaA.inp.'+str(it)
		print os.popen(cmd).read() 
   		cmd = 'cp Gf.out GfA.out.'+str(it)
		print os.popen(cmd).read() 
    		cmd = 'cp Sig.out SigA.out.'+str(it)
    		print os.popen(cmd).read() 
		cmd = 'cp ctqmc.log ctqmc.logA'+str(it)
                print os.popen(cmd).read() 

		# Constructing Delta.inp for B sublattice 
		params['mu'][0]=chemical_potential + delta
		CreateInputFileB(params)
		SigA = loadtxt('SigA.out.'+str(it))
    		cmd = "aprun -j 1 -n 24 -N 24 " +  "./hybri B " + str(chemical_potential) + " "+str(it)
                subprocess.call(cmd,shell=True,stdout=hyb_info,stderr=hyb_info)
   		hyb_info.flush()
		if it !=0: mix_deltaB(it,mix)
		

    		# Running ctqmc for B
    		print 'Running ---- qmc itt.B: ', it, '-----'
   		cmd = "aprun -j 1 -n 24 -N 24 "+params['exe'][0]+'  PARAMS > nohup_imp.out 2>&1 '
    		subprocess.call(cmd,shell=True,stdout=fh_info,stderr=fh_info)
   		fh_info.flush()
		
		#saving B
		cmd = 'cp Delta.inp DeltaB.inp.'+str(it)
		print os.popen(cmd).read() 
   		cmd = 'cp Gf.out GfB.out.'+str(it)
		print os.popen(cmd).read() # copying Gf
    		cmd = 'cp Sig.out SigB.out.'+str(it)
    		print os.popen(cmd).read() # copying Sig
		cmd = 'cp ctqmc.log ctqmc.logB'+str(it)
                print os.popen(cmd).read() # copying Sig
		
		#adjusting the chemical potential
		with open('ctqmc.logA'+str(it), 'r') as verFile:
    			for i, line in enumerate(verFile):
				if "nf" in line:
					text = line.split()
					text1 = text[0].split("=")
					nfA = float(text1[1])
				if "<<Sz>>" in line:
					text = line.split()
					text1 = text[3].split("=")
					SzA = float(text1[1])
		with open('ctqmc.logB'+str(it), 'r') as verFile:
    			for i, line in enumerate(verFile):
				if "nf" in line:
					text = line.split()
					text1 = text[0].split("=")
					nfB = float(text1[1])
				if "<<Sz>>" in line:
					text = line.split()
					text1 = text[3].split("=")
					SzB = float(text1[1])
		
		#if abs(nfA+nfB-2.0)>0.001 and l>=update_step and it<(Niter-update_end) and it> update_start:
		#	chemical_potential=chemical_potential-coeffi*(nfA+nfB-2.0)
		#	l=0
		print 'nfA,nfB,nfA+nfB',nfA,nfB,nfA+nfB
		print 'SzA,SzB,ms,mf',SzA,SzB,abs(SzA - SzB)/2.,abs(SzA + SzB)/2.
		print 'chemical potential',chemical_potential
	#all the iteration for given U is completed here
	cmd = 'cp *.out* *.dat* *.log* *.imp* PARAMS* *.cix* *.inp* *.00* %s/'%directory
   	print os.popen(cmd).read()

#TODO keep on appending the file that would be edited for hybf's output to check that everything is going write.
#TODO check how  subprocess work
