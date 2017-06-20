import pylab
import matplotlib.pyplot as plt
import numpy as np
import datetime


U_list  =[0.5,1.0,1.3] 
it_list =[79, 79, 79 ]
t = 0.5
beta = 100.0
for i,U in enumerate(U_list):
	print "U======",U
	data = np.loadtxt("GfA.out.%s.%s"%(it_list[i],U_list[i]))
	plt.plot(data[:,0], data[:,2],'-o',label = r"$\uparrow$")
	plt.plot(data[:,0], data[:,4],'-o',label = r"$\downarrow$")

	plt.xlim(0,10)
	plt.text(2.5,-0.5,r"$U=%st, \beta t = %s$"%(U/t, beta*t),size=20)
	plt.legend()
	plt.xlabel(r"$w_n$",size=20)
	plt.ylabel(r"$ImG_A$",size=20)
	plt.savefig("ImGA_U%s.eps"%U,format="eps")

	plt.show()

for i,U in enumerate(U_list):
	print "U======",U
	data = np.loadtxt("SigA.out.%s.%s"%(it_list[i],U_list[i]))
	plt.plot(data[:,0], data[:,2],'-o',label = r"$\uparrow$")
	plt.plot(data[:,0], data[:,4],'-o',label = r"$\downarrow$")

	plt.xlim(0,10)
	plt.text(2.5,-0.5,r"$U=%st, \beta t = %s$"%(U/t, beta*t),size=20)
	plt.legend()
	plt.xlabel(r"$w_n$",size=20)
	plt.ylabel(r"$Im\Sigma_A$",size=20)
	plt.savefig("ImSigA_U%s.eps"%U,format="eps")

	plt.show()

