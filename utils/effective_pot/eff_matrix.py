import numpy as np
import sys
import subprocess
from scipy.special import binom

ctcf_file = sys.argv[1]

command = "grep -c ^ "+(ctcf_file)
N=int(subprocess.check_output(command, shell=True))

with open(ctcf_file) as f:
    ctcf = [line.rstrip() for line in f]

pk=1
pl=1
Ueff = np.zeros((N,N))
for i in range(0,N-2):
    for j in range(i+2,N):
        sum_t=0
        for t in range(0, abs(i-j)-1+1):
            for k in range(0,t-1+1):
                if(int(ctcf[i+k])>0):
                    pk=0
                    break
                else:
                    pk=1
            for l in range(t+2, abs(i-j)+1):
                if(int(ctcf[i+l])<0):
                    pl=0
                    break
                else:
                    pl=1
            sum_t += binom(abs(i-j)-1, t)*pk*pl

        Ueff[i][j]=(1./2.57)*(-abs(i-j)*np.log(3e-13/(3e-13+1.28e-15))-np.log(2.28e-15/3e-13)-np.log(1./2.**(abs(i-j)-1)*sum_t))

Ueff_max=Ueff[np.isfinite(Ueff)].max()

for i in range(0,N-2):
   #print()
    for j in range(i+2,N):
        if(str(Ueff[i][j])=="inf"):
            print(i, j, 0)
        else:
            print(i, j, Ueff[i][j]-2.006362022611492)
