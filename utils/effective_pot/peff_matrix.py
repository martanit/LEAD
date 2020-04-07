import numpy as np
import sys
import subprocess
from scipy.special import binom

ctcf_file = sys.argv[1]

command = "grep -c ^ "+(ctcf_file)
N=int(subprocess.check_output(command, shell=True))

with open(ctcf_file) as f:
    ctcf = [line.rstrip() for line in f]

Ueff = np.zeros((N,N))
P_eff = np.zeros((N,N))
for i in range(0,N-2):
    for j in range(i+2,N):
        sum_t=0
        for t in range(0, abs(i-j)-1+1):
            pk=1
            pl=1
            for k in range(0,t-1+2):
                if(int(ctcf[i+k])>0):
                    pk=0
                    break
                else:
                    pk=1
            for l in range(t, abs(i-j)+1):
                if(int(ctcf[i+l])<0):
                    pl=0
                    break
                else:
                    pl=1
            sum_t += binom(abs(i-j)-1, t)*pk*pl

        P_eff[i][j] = (2.*2.28e-15/2e-13)*(2e-13/(2.*(2e-13+1.28e-15)))**(abs(i-j))*sum_t

for i in range(0,N-2):
    print()
    for j in range(i+2,N):
        print(i,j,P_eff[i][j])
