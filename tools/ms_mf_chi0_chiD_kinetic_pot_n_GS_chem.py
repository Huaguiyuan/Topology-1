import pylab
import matplotlib.pyplot as plt
import numpy as np
import datetime

beta = 100.0
t2 = 0.0
t1 = 0.5
delta = 0.0
last_it = [39,   39, 39,  39,  39,  39 ]
U_list =  [0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
EpotA = []
EpotB = []
EkinA =  []
EkinB = []
TRA = []
TRB = []
nfA,nfB = [],[]
chiS0A,chiS0B = [],[]
chiSDA,chiSDB = [],[]
szA,szB = [],[]
for l,U in enumerate(U_list):
	
	with open('ctqmc_%s.logA%s'%(U,last_it[l]), 'r') as verFile:
    		for i, line in enumerate(verFile):
			
	
			if "Epot" in line:
				text = line.split()
				text1 = text[0].split("=")
				EpotA.append(float(text1[1]))
				#print U,float(text1[1])
				text1 = text[1].split("=")
				EkinA.append(float(text1[1]))
				text1 = text[2].split("=")
				TRA.append(float(text1[1]))
			if "nf" in line:
				text = line.split()
				text1 = text[0].split("=")
				nfA.append(float(text1[1]))
				#print U,float(text1[1])
				text1 = text[1].split("=")
				chiS0A.append(float(text1[1]))
				text1 = text[2].split("=")
				chiSDA.append(float(text1[1]))
				text1 = text[3].split("=")
				szA.append(float(text1[1]))
	with open('ctqmc_%s.logB%s'%(U,last_it[l]), 'r') as verFile:
    		for i, line in enumerate(verFile):
			
			if "Epot" in line:
				text = line.split()
				text1 = text[0].split("=")
				EpotB.append(float(text1[1]))
				#print U,float(text1[1])
				text1 = text[1].split("=")
				EkinB.append(float(text1[1]))
				text1 = text[2].split("=")
				TRB.append(float(text1[1]))
			if "nf" in line:
				text = line.split()
				text1 = text[0].split("=")
				nfB.append(float(text1[1]))
				#print U,float(text1[1])
				text1 = text[1].split("=")
				chiS0B.append(float(text1[1]))
				text1 = text[2].split("=")
				chiSDB.append(float(text1[1]))
				text1 = text[3].split("=")
				szB.append(float(text1[1]))
				
				
			'''
			if "Epot" in line:
				print line[6:15]
			if "Ekin" in line:
				print line[23:32]
			if 'Tr(Sigma*G)' in line:
				print line[46:55]
			if 'nf' in line:
				print line[3:10]
			if 'chiS0' in line:
				print line[17:25]
			if 'chiD' in line:
				print line[31:39]
			if '<<Sz>>' in line:
				print line[47:59]
			'''
	#lA = pylab.loadtxt('ctqmc%s.logA'%U)
	#lB = pylab.loadtxt('ctqmc%s.logB'%U)
	#print lA 
#len(U_list),len(EpotA)	

#plt.plot(U_list,EpotA,label = 'Epot')
#plt.plot(U_list,EkinA,label = 'Ekin')
#plt.plot(U_list,TRA,label = r'$TrG\Sigma$')
#plt.plot(U_list,chiS0A,label = r'$\chi_0$')
#plt.plot(U_list,chiSDA,label = r'$\chi_D$')
n = (np.array(nfB) + np.array(nfA))/2.0
ms = abs(np.array(szB) - np.array(szA))/2.0
mf = abs(np.array(szB) + np.array(szA))/2.0
print "n",n
f1 = open('result1_t2%s_%s_%s_date%s.dat'%((t2/t1),delta/t1,beta,datetime.date.today()),'w')
print >> f1, '# U/t	nA	nB	(nA+nB)/2	ms=abs(szA-szB)/2.0	mf=abs(szA+szB)/2.0'
for i in range(len(U_list)):
	print U_list[i]
	print >>f1, U_list[i]/t1,		nfA[i],	nfB[i],	n[i],	ms[i],	mf[i] 
f1.close()

f1 = open('result2_t2%s_%s_%s_date%s.dat'%((t2/t1),delta/t1,beta,datetime.date.today()),'w')
print >>f1, '# U/t	EkinA	EkinB	EpotA	EpotB	chiS0A	chiS0B	chiSDA	chiSDB'
for i in range(len(U_list)):
	print >>f1, U_list[i]/t1,	EkinA[i],	EkinB[i],  EpotA[i],  EpotB[i],	chiS0A[i], chiS0B[i], chiSDA[i], chiSDB[i]
f1.close()


# plotting staggered magnetization.
plt.plot(np.array(U_list)/t1,abs(np.array(szA)- np.array(szB)),'--o',label = r'$m_s$')
plt.ylabel(r'$m_s$',size=20)
plt.tick_params(axis='y', labelsize=15)
plt.xlabel(r'$U/t$',size=20)
plt.tick_params(axis='y', labelsize=15)
plt.text(1,0.2,r'$\Delta/t=%s,\beta=%s,t^\prime=%st$'%(delta/t1,beta,t2/t1),size=20)
plt.savefig("IHM_beta%s_t_prime_%s_Mf_vs_U_half_filling_date%s.pdf"%(beta,(t2/t1),datetime.date.today()))
plt.legend()
#plt.xlim(0,10)
plt.show()

# plotting ferro magnetization.
plt.plot(np.array(U_list)/t1,abs(np.array(szA)+np.array(szB)),'--o',label = r'$m_f$')
plt.ylabel(r'$m_f$',size=20)
plt.tick_params(axis='y', labelsize=15)
plt.xlabel(r'$U/t$',size=20)
plt.tick_params(axis='y', labelsize=15)
plt.text(1,0.2,r'$\Delta/t=%s,\beta=%s,t^\prime=%st$'%(delta/t1,beta,t2/t1),size=20)
plt.savefig("IHM_beta%s_t_prime_%s_Ms_vs_U_half_filling_date%s.pdf"%(beta,(t2/t1),datetime.date.today()))
plt.legend()
#plt.xlim(0,10)
plt.show()

# plotting occupency.
plt.plot(np.array(U_list)/t1,np.array(nfA),'--o',label = r'$n_{nfA}$')
plt.plot(np.array(U_list)/t1,np.array(nfB),'--o',label = r'$n_{nfB}$')
plt.plot(np.array(U_list)/t1,(np.array(nfB) + np.array(nfA))/2.0,'--o',label = r'$n_{ave}$')
plt.plot(np.array(U_list)/t1,abs(np.array(nfB) - np.array(nfA))/2.0,'--o',label = r'$\delta n$')
#plt.ylabel(r'$m_s$',size=20)
plt.tick_params(axis='y', labelsize=15)
plt.xlabel(r'$U/t$',size=20)
plt.tick_params(axis='y', labelsize=15)
plt.text(1,0.2,r'$\Delta/t=%s,\beta=%s,t^\prime=%st$'%(delta/t1,beta,t2/t1),size=20)
plt.savefig("IHM_beta%s_t_prime_%s_n_vs_U_half_filling_date%s.pdf"%(beta,(t2/t1),datetime.date.today()))
plt.legend()
#plt.xlim(0,10)
plt.show()

# plotting kinetic energy.
plt.plot(np.array(U_list)/t1,np.array(EkinA),'--o',label = r'$E_{kA}$')
plt.plot(np.array(U_list)/t1,np.array(EkinB),'--o',label = r'$E_{kB}$')
plt.plot(np.array(U_list)/t1,(np.array(EkinB)+np.array(EkinA))/2.0,'--o',label = r'$E_{kave}$')
plt.ylabel(r'$E_k$',size=20)
plt.tick_params(axis='y', labelsize=15)
plt.xlabel(r'$U/t$',size=20)
plt.tick_params(axis='y', labelsize=15)
plt.text(1,-0.34,r'$\Delta/t=%s,\beta=%s,t^\prime=%st$'%(delta/t1,beta,t2/t1),size=20)
plt.savefig("IHM_beta%s_t_prime_%s_Ekin_vs_U_half_filling_date%s.pdf"%(beta,(t2/t1),datetime.date.today()))
plt.legend()
plt.show()
# plotting potential energy.
plt.plot(np.array(U_list)/t1,np.array(EpotA),'--o',label = r'$E_{potA}$')
plt.plot(np.array(U_list)/t1,np.array(EpotB),'--o',label = r'$E_{potB}$')
plt.plot(np.array(U_list)/t1,(np.array(EpotB)+np.array(EpotA))/2.0,'--o',label = r'$E_{potave}$')
plt.ylabel(r'$E_{pot}$',size=20)
plt.tick_params(axis='y', labelsize=15)
plt.xlabel(r'$U/t$',size=20)
plt.tick_params(axis='y', labelsize=15)
plt.text(1,0.2,r'$\Delta/t=%s,\beta=%s,t^\prime=%st$'%(delta/t1,beta,t2/t1),size=20)
plt.savefig("IHM_beta%s_t_prime_%s_Epot_vs_U_half_filling_date%s.pdf"%(beta,(t2/t1),datetime.date.today()))
plt.legend()
#plt.xlim(0,10)
plt.show()

# plotting sciS0.
plt.plot(np.array(U_list)/t1,np.array(chiS0A),'--o',label = r'$\chi_{0A}$')
plt.plot(np.array(U_list)/t1,np.array(chiS0B),'--o',label = r'$\chi_{0B}$')
plt.plot(np.array(U_list)/t1,(np.array(chiS0A)+np.array(chiS0B))/2.0,'--o',label = r'$\chi_{0}$')
plt.ylabel(r'$\chi_{0}$',size=20)
plt.tick_params(axis='y', labelsize=15)
plt.xlabel(r'$U/t$',size=20)
plt.tick_params(axis='y', labelsize=15)
plt.text(1,6,r'$\Delta/t=%s,\beta=%s,t^\prime=%st$'%(delta/t1,beta,t2/t1),size=20)
plt.savefig("IHM_beta%s_t_prime_%s_chiS0_vs_U_half_filling_date%s.pdf"%(beta,(t2/t1),datetime.date.today()))
plt.legend()
#plt.xlim(0,10)
plt.show()

# plotting sciSD.
plt.plot(np.array(U_list)/t1,np.array(chiSDA),'--o',label = r'$\chi_{DA}$')
plt.plot(np.array(U_list)/t1,np.array(chiSDB),'--o',label = r'$\chi_{DB}$')
plt.plot(np.array(U_list)/t1,(np.array(chiSDA)+np.array(chiSDB))/2.0,'--o',label = r'$\chi_{D}$')
plt.ylabel(r'$\chi_{D}$',size=20)
plt.tick_params(axis='y', labelsize=15)
plt.xlabel(r'$U/t$',size=20)
plt.tick_params(axis='y', labelsize=15)
plt.text(1,0.1,r'$\Delta/t=%s,\beta=%s,t^\prime=%st$'%(delta/t1,beta,t2/t1),size=20)
plt.savefig("IHM_beta%s_t_prime_%s_chiSD_vs_U_half_filling_date%s.pdf"%(beta,(t2/t1),datetime.date.today()))
plt.legend()
#plt.xlim(0,10)
plt.show()
