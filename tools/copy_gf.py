#from scipy import *
import numpy as np
import os
import fnmatch




U = 2.1
N=79
import os
for root, dirs, files in os.walk(".", topdown=False):
    for name in files:
        #print(os.path.join(root, name))
	#print "name:",name
    	#if (name=="initials*"):
	#if (fnmatch.fnmatch(name, "GfA.out.%s*"%N) or fnmatch.fnmatch(name, "GfB.out.%s*"%N) or fnmatch.fnmatch(name, "SigB.out.%s*"%N) or fnmatch.fnmatch(name, "SigB.out.%s*"%N) or fnmatch.fnmatch(name, "DeltaA.inp.%s*"%N) or fnmatch.fnmatch(name, "DeltaB.inp.%s*"%N)):
        #print name, ["GfA.out.%s"%N,"GfB.out.%s"%N, "SigB.out.%s"%N, "SigA.out.%s"%N, "DeltaA.inp.%s"%N, "DeltaB.inp.%s"%N]
	if name in ["GfA.out.%s"%N,"GfB.out.%s"%N, "SigB.out.%s"%N, "SigA.out.%s"%N, "DeltaA.inp.%s"%N, "DeltaB.inp.%s"%N]:
		name1 = name +".%s"%U
    		cmd = "cp ./%s ./%s"%(name,name1)
                print name1
    		print os.popen(cmd).read()
		#print os.path.splitext(name)[0]+"_%s"%U + ".dat"
    		cmd = "cp %s  ../analysis/green"%name1
		print os.popen(cmd).read()
	if name in ["ctqmc.logA%s"%N,"ctqmc.logB%s"%N]:
		#3name1 = name +".%s"%U
		name1 = os.path.splitext(name)[0]+"_%s"%U + os.path.splitext(name)[1]
    		cmd = "cp ./%s ./%s"%(name,name1)
                print name1
    		print os.popen(cmd).read()
		#print os.path.splitext(name)[0]+"_%s"%U + ".dat"
    		cmd = "cp %s  ../analysis"%name1
		print os.popen(cmd).read()
	if name == "nohup_imp.out.000":
		name1 = name +"_%s"%U
		cmd = "cp ./%s ./%s"%(name,name1)
                print name1
    		print os.popen(cmd).read()
		cmd = "cp %s  ../analysis"%name1
		print os.popen(cmd).read()


