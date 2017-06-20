import numpy as np

beta=100.0
N = 1025
SigA = open("SigA.out.0",'w')
SigB = open("SigB.out.0",'w')
for i in range(N):
        w = (2*i-1)*np.pi/beta
        print >> SigA, w, 0.0, 0.0, 0.0, 0.0
        print >> SigB, w, 0.0, 0.0, 0.0, 0.0
