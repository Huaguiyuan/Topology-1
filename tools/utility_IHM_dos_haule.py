
import numpy as np
from numpy import *
#import matplotlib.pyplot as plt

class green:
	
	def __init__(self,beta,delta):
		self.beta = beta
		self.N = 1025
		self.delta = delta
		self.G0Au = np.zeros(self.N,complex)	
		self.G0Bu = np.zeros(self.N,complex)
		self.G0Ad = np.zeros(self.N,complex)	
		self.G0Bd = np.zeros(self.N,complex)
		self.self_energyAu = np.zeros(self.N,complex)	
		self.self_energyBu = np.zeros(self.N,complex)
		self.self_energyAd = np.zeros(self.N,complex)	
		self.self_energyBd = np.zeros(self.N,complex)
		self.wn = np.zeros(self.N)
		self.GAu = np.zeros(self.N,complex)
		self.GBu = np.zeros(self.N,complex)
		self.GAd = np.zeros(self.N,complex)
		self.GBd = np.zeros(self.N,complex)
		self.DAu = np.zeros(self.N,complex)
		self.DBu = np.zeros(self.N,complex)
		self.DAd = np.zeros(self.N,complex)
		self.DBd = np.zeros(self.N,complex)
		for i in range(self.N):
			self.wn[i] = (2*i + 1) * np.pi/self.beta
				
		
	def non_iteracting_G_IHM(self,dos_file, mu):
		E = loadtxt(dos_file)
		length = len(E)
		step = abs(E[1][0] - E[0][0])
		#W=np.arange(-1.99,2.001,de)
		for k in range(self.N):
			alphaA = 1j*self.wn[k] - self.delta +mu
			alphaB = 1j*self.wn[k] + self.delta +mu
			gwA = 0.0
			gwB = 0.0
			for i in range(length):
				l = 0.0
				m = E[i][0]
				rho= E[i][1]
				if i ==0 or i == length: 
					gwA = gwA + rho*(alphaB-l)/((alphaA-l)*(alphaB-l)-m**2)
					gwB = gwB + rho*(alphaA-l)/((alphaA-l)*(alphaB-l)-m**2)
				elif i%2 ==0 : 
					gwA = gwA + 2*rho*(alphaB-l)/((alphaA-l)*(alphaB-l)-m**2)
					gwB = gwB + 2*rho*(alphaA-l)/((alphaA-l)*(alphaB-l)-m**2)
				else : 
					gwA = gwA + 4*rho*(alphaB-l)/((alphaA-l)*(alphaB-l)-m**2)
					gwB = gwB + 4*rho*(alphaA-l)/((alphaA-l)*(alphaB-l)-m**2)
				
			self.GAu[k] = gwA * (step/3.)
			self.GBu[k] = gwB * (step/3.)
			self.GAd[k] = gwA * (step/3.)
			self.GBd[k] = gwB * (step/3.)
	def hybri_function_IHM_dosA(self,fDelta,dos_file,mu):
		E = loadtxt(dos_file)
		length = len(E)
		step = abs(E[1][0] - E[0][0])
		f = open(fDelta,'w')
		#W=np.arange(-1.99,2.001,de)
		for k in range(self.N):
			alphaAu = 1j*self.wn[k] - self.delta +mu - self.self_energyAu[k]
			alphaBu = 1j*self.wn[k] + self.delta +mu - self.self_energyBu[k]
			gwAu = 0.0
			#gwBu = 0.0
			alphaAd = 1j*self.wn[k] - self.delta +mu - self.self_energyAd[k]
			alphaBd = 1j*self.wn[k] + self.delta +mu - self.self_energyBd[k]
			gwAd = 0.0
			#gwBd = 0.0
			for i in range(length):
				l = 0.0
				m = E[i][0]
				rho= E[i][1]
				if i ==0 or i == length: 
					gwAu = gwAu + rho*(alphaBu-l)/((alphaAu-l)*(alphaBu-l)-m**2)
					#gwBu = gwBu + rho*(alphaAu-l)/((alphaAu-l)*(alphaBu-l)-m**2)
					gwAd = gwAd + rho*(alphaBd-l)/((alphaAd-l)*(alphaBd-l)-m**2)
					#gwBd = gwBd + rho*(alphaAd-l)/((alphaAd-l)*(alphaBd-l)-m**2)
				elif i%2 ==0 : 
					gwAu = gwAu + 2*rho*(alphaBu-l)/((alphaAu-l)*(alphaBu-l)-m**2)
					#gwBu = gwBu + 2*rho*(alphaAu-l)/((alphaAu-l)*(alphaBu-l)-m**2)
					gwAd = gwAd + 2*rho*(alphaBd-l)/((alphaAd-l)*(alphaBd-l)-m**2)
					#gwBd = gwBd + 2*rho*(alphaAd-l)/((alphaAd-l)*(alphaBd-l)-m**2)
				else : 
					gwAu = gwAu + 4*rho*(alphaBu-l)/((alphaAu-l)*(alphaBu-l)-m**2)
					#gwBu = gwBu + 2*rho*(alphaAu-l)/((alphaAu-l)*(alphaBu-l)-m**2)
					gwAd = gwAd + 4*rho*(alphaBd-l)/((alphaAd-l)*(alphaBd-l)-m**2)
					#gwBd = gwBd + 2*rho*(alphaAd-l)/((alphaAd-l)*(alphaBd-l)-m**2)
				
			self.GAu[k] = gwAu * (step/3.)
			#self.GBu[k] = gwBu * (step/3.)
			self.GAd[k] = gwAd * (step/3.)
			#self.GBd[k] = gwBd * (step/3.)
			#for k in range(self.N):	
			self.G0Au[k] = 1/(1/self.GAu[k] + self.self_energyAu[k])
			self.G0Ad[k] = 1/(1/self.GAd[k] + self.self_energyAd[k])
			self.DAu[k] = 1j*self.wn[k] + mu - self.delta - 1./self.G0Au[k]
			self.DAd[k] = 1j*self.wn[k] + mu - self.delta - 1./self.G0Ad[k]
			print >> f, self.wn[k], self.DAu[k].real, self.DAu[k].imag, self.DAd[k].real, self.DAd[k].imag
		f.close()

	def hybri_function_IHM_dosB(self,fDelta,dos_file,mu):
		E = loadtxt(dos_file)
		length = len(E)
		step = abs(E[1][0] - E[0][0])
		f = open(fDelta, 'w')
		#W=np.arange(-1.99,2.001,de)
		for k in range(self.N):
			alphaAu = 1j*self.wn[k] - self.delta +mu - self.self_energyAu[k]
			alphaBu = 1j*self.wn[k] + self.delta +mu	- self.self_energyBu[k]
			#gwAu = 0.0
			gwBu = 0.0
			alphaAd = 1j*self.wn[k] - self.delta +mu - self.self_energyAd[k]
			alphaBd = 1j*self.wn[k] + self.delta +mu	- self.self_energyBd[k]
			#gwAd = 0.0
			gwBd = 0.0
			for i in range(length):
				l = 0.0
				m = E[i][0]
				rho= E[i][1]
				if i ==0 or i == length: 
					#gwAu = gwAu + rho*(alphaBu-l)/((alphaAu-l)*(alphaBu-l)-m**2)
					gwBu = gwBu + rho*(alphaAu-l)/((alphaAu-l)*(alphaBu-l)-m**2)
					#gwAd = gwAd + rho*(alphaBd-l)/((alphaAd-l)*(alphaBd-l)-m**2)
					gwBd = gwBd + rho*(alphaAd-l)/((alphaAd-l)*(alphaBd-l)-m**2)
				elif i%2 ==0 : 
					#gwAu = gwAu + 2*rho*(alphaBu-l)/((alphaAu-l)*(alphaBu-l)-m**2)
					gwBu = gwBu + 2*rho*(alphaAu-l)/((alphaAu-l)*(alphaBu-l)-m**2)
					#gwAd = gwAd + 2*rho*(alphaBd-l)/((alphaAd-l)*(alphaBd-l)-m**2)
					gwBd = gwBd + 2*rho*(alphaAd-l)/((alphaAd-l)*(alphaBd-l)-m**2)
				else : 
					#gwAu = gwAu + 2*rho*(alphaBu-l)/((alphaAu-l)*(alphaBu-l)-m**2)
					gwBu = gwBu + 4*rho*(alphaAu-l)/((alphaAu-l)*(alphaBu-l)-m**2)
					#gwAd = gwAd + 2*rho*(alphaBd-l)/((alphaAd-l)*(alphaBd-l)-m**2)
					gwBd = gwBd + 4*rho*(alphaAd-l)/((alphaAd-l)*(alphaBd-l)-m**2)
				
			#self.GAu[k] = gwAu * (step/3.)
			self.GBu[k] = gwBu * (step/3.)
			#self.GAd[k] = gwAd * (step/3.)
			self.GBd[k] = gwBd * (step/3.)
		#for k in range(self.N):	
			self.G0Bu[k] = 1/(1/self.GBu[k] + self.self_energyBu[k])
			self.G0Bd[k] = 1/(1/self.GBd[k] + self.self_energyBd[k])
			self.DBu[k] = 1j*self.wn[k] + mu + self.delta - 1./self.G0Bu[k]
			self.DBd[k] = 1j*self.wn[k] + mu + self.delta - 1./self.G0Bd[k]
			
			print >> f,self.wn[k], self.DBu[k].real, self.DBu[k].imag, self.DBd[k].real, self.DBd[k].imag
		f.close()
		
	def mix_deltaA(self,ita,mix):
		dold = loadtxt('DeltaA.inp.'+str(ita-1))
		dnew = loadtxt('Delta.inp')
		delta = open('Delta.inp','w')
		data = mix*dnew + (1-mix)*dold
		for i in range(self.N):
			print >> delta, dold[i][0], data[i][1], data[i][2], data[i][3], data[i][4]
		delta.close()
	
	def mix_deltaB(self,ita,mix):
		dold = loadtxt('DeltaB.inp.'+str(ita-1))
		dnew = loadtxt('Delta.inp')
		delta = open('Delta.inp','w')
		data = mix*dnew + (1-mix)*dold
		for i in range(self.N):
			print >> delta, dold[i][0], data[i][1], data[i][2], data[i][3], data[i][4]
		delta.close()


	def CreateInputFileA(self,params):
	    " Creates input file (PARAMS) for CT-QMC solver"
	    fA = open('PARAMSA', 'w')
	    f = open('PARAMS', 'w')
	    print >> fA, '# Input file for continuous time quantum Monte Carlo'
	    print >> f, '# Input file for continuous time quantum Monte Carlo'
	    for p in params:
		print >> f, p, params[p][0], '\t', params[p][1]
		print >> fA, p, params[p][0], '\t', params[p][1]
	    f.close()

	def CreateInputFileB(self,params):
	    " Creates input file (PARAMS) for CT-QMC solver"
	    f = open('PARAMS', 'w')
	    fB = open('PARAMSB', 'w')
	    print >> f, '# Input file for continuous time quantum Monte Carlo'
	    print >> fB, '# Input file for continuous time quantum Monte Carlo'
	    for p in params:
		print >> f, p, params[p][0], '\t', params[p][1]
		print >> fB, p, params[p][0], '\t', params[p][1]
	    f.close()

	def create_cix_input(self,params,icix):
		f = open(params['cix'][0], 'w')
		print >> f, icix
		f.close()
		
	def non_iteracting_G_IHM_bethett2(self, mu,t,t2):
		de = 0.01
		E=np.arange(-1.9999,1.9999,de)
		eta=0.001
		length = len(E)
		#W=np.arange(-1.99,2.001,de)
		for k in range(self.N):
			alphaA = 1j*self.wn[k] - self.delta +mu
			alphaB = 1j*self.wn[k] + self.delta +mu
			gwA = 0.0
			gwB = 0.0
			for i,theta in enumerate(E):
				l = t2*(theta**2 -1)
				m = (theta*t)
				rho= np.sqrt(4-theta**2)/(2*np.pi)
				if i ==0 or i == length: 
					gwA = gwA + rho*(alphaB-l)/((alphaA-l)*(alphaB-l)-m**2)
					gwB = gwB + rho*(alphaA-l)/((alphaA-l)*(alphaB-l)-m**2)
				elif i%2 ==0 : 
					gwA = gwA + 2*rho*(alphaB-l)/((alphaA-l)*(alphaB-l)-m**2)
					gwB = gwB + 2*rho*(alphaA-l)/((alphaA-l)*(alphaB-l)-m**2)
				else : 
					gwA = gwA + 4*rho*(alphaB-l)/((alphaA-l)*(alphaB-l)-m**2)
					gwB = gwB + 4*rho*(alphaA-l)/((alphaA-l)*(alphaB-l)-m**2)
				
			self.G0A[k] = gwA * (de/3.)
			self.G0B[k] = gwB * (de/3.)
	
	def iteracting_lattice_G_IHM(self,dos_file, mu):
		E = loadtxt(dos_file)
		length = len(E)
		step = abs(E[1][0] - E[0][0])
		#W=np.arange(-1.99,2.001,de)
		for k in range(self.N):
			alphaAu = 1j*self.wn[k] - self.delta +mu - self.self_energyAu[k]
			alphaBu = 1j*self.wn[k] + self.delta +mu	- self.self_energyBu[k]
			gwAu = 0.0
			gwBu = 0.0
			alphaAd = 1j*self.wn[k] - self.delta +mu - self.self_energyAd[k]
			alphaBd = 1j*self.wn[k] + self.delta +mu - self.self_energyBd[k]
			gwAd = 0.0
			gwBd = 0.0
			for i in range(length):
				l = 0.0
				m = E[i][0]
				rho= E[i][1]
				if i ==0 or i == length: 
					gwAu = gwAu + rho*(alphaBu-l)/((alphaAu-l)*(alphaBu-l)-m**2)
					gwBu = gwBu + rho*(alphaAu-l)/((alphaAu-l)*(alphaBu-l)-m**2)
					gwAd = gwAd + rho*(alphaBd-l)/((alphaAd-l)*(alphaBd-l)-m**2)
					gwBd = gwBd + rho*(alphaAd-l)/((alphaAd-l)*(alphaBd-l)-m**2)
				elif i%2 ==0 : 
					gwAu = gwAu + 2*rho*(alphaBu-l)/((alphaAu-l)*(alphaBu-l)-m**2)
					gwBu = gwBu + 2*rho*(alphaAu-l)/((alphaAu-l)*(alphaBu-l)-m**2)
					gwAd = gwAd + 2*rho*(alphaBd-l)/((alphaAd-l)*(alphaBd-l)-m**2)
					gwBd = gwBd + 2*rho*(alphaAd-l)/((alphaAd-l)*(alphaBd-l)-m**2)
				else : 
					gwAu = gwAu + 2*rho*(alphaBu-l)/((alphaAu-l)*(alphaBu-l)-m**2)
					gwBu = gwBu + 2*rho*(alphaAu-l)/((alphaAu-l)*(alphaBu-l)-m**2)
					gwAd = gwAd + 2*rho*(alphaBd-l)/((alphaAd-l)*(alphaBd-l)-m**2)
					gwBd = gwBd + 2*rho*(alphaAd-l)/((alphaAd-l)*(alphaBd-l)-m**2)
				
			self.GAu[k] = gwAu * (step/3.)
			self.GBu[k] = gwBu * (step/3.)
			self.GAd[k] = gwAd * (step/3.)
			self.GBd[k] = gwBd * (step/3.)
	
		
	
	# NON interacting lattice green's function for sigle site problem for the time being A site was chosen.	
	def non_interating_G(self,dos_file,mu):
		E = loadtxt(dos_file)
		length = len(E)
		step = abs(E[1][0] - E[0][0])
		for k in range(self.N):
			gw = 0.0
			for i in range(length):
				if i ==0 or i == length: gw = gw + E[i][1]/ (1j*self.wn[k] - mu- E[i][0])
				elif i%2 ==0 : gw = gw + (E[i][1] * 2)/ (1j*self.wn[k] - mu- E[i][0])
				else : gw = gw + (E[i][1] * 4)/ (1j*self.wn[k] - mu- E[i][0])
			self.G0Au[k] = gw * (step/3.)

	# interacting lattice green's function for sigle site problem for the time being A site was chosen.	
	def interating_lattice_G(self,dos_file,mu):
		
		E = loadtxt(dos_file)
		length = len(E)
		step = abs(E[1][0] - E[0][0])
		for k in range(self.N):
			gw = 0.0
			for i in range(length):
				if i ==0 or i == length: gw = gw + E[i][1]/ (1j*self.wn[k] - mu- E[i][0]- self.self_energyAu[k])
				elif i%2 ==0 : gw = gw + (E[i][1] * 2)/ (1j*self.wn[k] - mu- E[i][0]- self.self_energyAu[k])
				else : gw = gw + (E[i][1] * 4)/ (1j*self.wn[k] - mu- E[i][0]- self.self_energyAu[k])
			self.GAu[k] = gw * (step/3.)

	# wiess function for sigle site problem for the time being A site was chosen.
	def wies_field(self,):
		for k in range(self.N):	
			self.G0Au[k] = 1/(1/self.GAu[k] + self.self_energyAu[k])
			self.DAu[k] = 1j*self.wn[k] - 1./self.GAu[k]
			self.DAd[k] = 1j*self.wn[k] - 1./self.G0Au[k]
		

