import pylab
import matplotlib.pyplot as plt
import numpy as np
import datetime


U_list  =[0.5] 
it_list =[79, 79, 79 ]
t = 0.5
beta = 100.0
for i,U in enumerate(U_list):
	print "U======",U
	#data = np.loadtxt("GfA.out.%s.%s"%(it_list[i],U_list[i]))
        data = np.loadtxt("Gf.out")
	plt.plot(data[:,0], data[:,2],'-o',label = r"$\uparrow$")
	plt.plot(data[:,0], data[:,8],'-o',label = r"$\downarrow$")
        plt.plot(data[:,0], data[:,4],'-o',label = r"$\uparrow \downarrow$")
	plt.plot(data[:,0], data[:,6],'-o',label = r"$\downarrow \uparrow$")

	plt.xlim(0,10)
	plt.text(2.5,-0.5,r"$U=%st, \beta t = %s$"%(U/t, beta*t),size=20)
	plt.legend()
	plt.xlabel(r"$w_n$",size=20)
	plt.ylabel(r"$ImG$",size=20)
	plt.savefig("ImG_U%s.eps"%U,format="eps")

	plt.show()

for i,U in enumerate(U_list):
	print "U======",U
	#data = np.loadtxt("SigA.out.%s.%s"%(it_list[i],U_list[i]))
        data = np.loadtxt("Sig.out")
	plt.plot(data[:,0], data[:,2],'-o',label = r"$\uparrow$")
	plt.plot(data[:,0], data[:,8],'-o',label = r"$\downarrow$")
        plt.plot(data[:,0], data[:,4],'-o',label = r"$\uparrow \downarrow$")
	plt.plot(data[:,0], data[:,6],'-o',label = r"$\downarrow \uparrow$")

	plt.xlim(0,10)
	plt.text(2.5,-0.5,r"$U=%st, \beta t = %s$"%(U/t, beta*t),size=20)
	plt.legend()
	plt.xlabel(r"$w_n$",size=20)
	plt.ylabel(r"$Im\Sigma$",size=20)
	plt.savefig("ImSigU%s.eps"%U,format="eps")

	plt.show()

