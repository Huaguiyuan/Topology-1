import matplotlib.pyplot as plt
import numpy as np
import datetime
from scipy import *
#U_list=[0.0,0.2,0.4,0.5,0.7,1.0,1.3,1.5,1.625,1.75,1.8125,1.875,2.0,2.1,2.2,2.3,2.4,2.5,3.0,3.5,4.0,5.0,6.0,7.0,7.5,8.0,8.5,9.0,9.5]
it_list = [49,  49, 69,   79, 59,  49,  49,   59, 119,  59]
U_list  = [0.0,0.5, 0.6, 0.8, 0.9, 1.0, 1.2]#,0.5,1.0,1.5,1.625,1.6875,1.75,2.0,2.5]
delta=0.0
beta=50
t=0.5
t2=0.0*t
g0tA_up=[]
g0tB_up=[]
g0tA_down=[]
g0tB_down=[]
for i,U in enumerate(U_list):

	GA = loadtxt('GfA.out.%s.%s'%(it_list[i],U))
        GB = loadtxt('GfB.out.%s.%s'%(it_list[i],U))
        Au = 0.0
        Ad = 0.0
        Bu = 0.0
        Bd = 0.0

        for i in range(len(GA[:,0])):
		Au = Au + 2*( GA[i,1]*np.cos(GA[i,0]*0.5*beta) - GA[i,2]*np.sin(GA[i,0]*0.5*beta) )
                Ad = Ad + 2*( GA[i,3]*np.cos(GA[i,0]*0.5*beta) - GA[i,4]*np.sin(GA[i,0]*0.5*beta) )
                Bu = Bu + 2*( GB[i,1]*np.cos(GB[i,0]*0.5*beta) - GB[i,2]*np.sin(GB[i,0]*0.5*beta) )
                Bd = Bd + 2*( GB[i,3]*np.cos(GB[i,0]*0.5*beta) - GB[i,4]*np.sin(GB[i,0]*0.5*beta) )
	g0tA_up.append(Au)
        g0tA_down.append(Ad)
	g0tB_up.append(Bu)
        g0tB_down.append(Bd)
'''
	
	glist = [ GfImTime(indices = al, beta = beta) for a,al in gf_struct]
	g0tA = BlockGf(name_list = a_list, block_list = glist, make_copies=True, name="g0t")
	g0tB = BlockGf(name_list = a_list, block_list = glist, make_copies=True, name="g0t")
	g0tA['up'] <<= InverseFourier(gA['up'])
	g0tB['up'] <<= InverseFourier(gB['up'])
	g0tA['down'] <<= InverseFourier(gA['down'])
	g0tB['down'] <<= InverseFourier(gB['down'])
	#g0t['down'] <<= InverseFourier(g['down'])
	g0tA_up.append(g0tA['up'].data[5000][0][0])
	g0tB_up.append(g0tB['up'].data[5000][0][0])
	g0tA_down.append(g0tA['down'].data[5000][0][0])
	g0tB_down.append(g0tB['down'].data[5000][0][0])
R=HDFArchive("IHM_Uall"+"_beta%s"%beta + "_delta%s"%delta + "_mu%s"%0.0+"G_beta_by_2.h5")
R['g0tA_up']=g0tA_up
R['g0tB_up']=g0tB_up
R['g0tA_down']=g0tA_down
R['g0tB_down']=g0tB_down
t=0.5
'''
#beta = 1
plt.plot(np.array(U_list)/t,np.array(g0tA_up),'--o',label=r'$\tilde{A}_{A\uparrow}$',markersize=6,linewidth=2.5)
plt.plot(np.array(U_list)/t,np.array(g0tB_up),'--*',label=r'$\tilde{A}_{B\uparrow}$',markersize=6,linewidth=2.5)
plt.plot(np.array(U_list)/t,np.array(g0tA_down),'--s',label=r'$\tilde{A}_{A\downarrow}$',markersize=6,linewidth=2.5)
plt.plot(np.array(U_list)/t,np.array(g0tB_down),'--^',label=r'$\tilde{A}_{B\downarrow}$',markersize=6,linewidth=2.5)
plt.plot([0.775/t,0.775/t],[0.0,0.25],'--k')
plt.legend(prop={'size':27},frameon=False,labelspacing=0.01,borderaxespad=0.01,columnspacing=0.01,handletextpad=0.005)
plt.gcf().subplots_adjust(bottom=0.17)
plt.gcf().subplots_adjust(left=0.20)
plt.ylabel(r'$\tilde{A}(w=0)$',size=30)
plt.tick_params(axis='y', labelsize=22)
plt.xlabel(r'$U/t$',size=30)
plt.tick_params(axis='x', labelsize=22)
plt.text(0.5,0.1,r'$\Delta/t=%s,t2=%st$'%(delta/t,t2/t),size=25)
plt.savefig("IHM_G_beta_by_2_ms_vs_beta%s_delta%s_half_filing_bate%s.eps"%(beta,delta,datetime.date.today()))
plt.show()

f = open('A_0_%s_%s_%s.dat'%(delta/t,beta,t2/t),'w')
print >> f, "#U/t A_up  A_dn B_up  B_dn"
print >> f, "#Uc/t = 1.55"
for i in range(len(U_list)):
	print >> f, U_list[i]/t , g0tA_up[i], g0tA_down[i], g0tB_up[i], g0tB_down[i]
f.close()


	

