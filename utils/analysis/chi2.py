import sys
import numpy as np

exp_in=open(sys.argv[1], 'r')
sim_in=open(sys.argv[2], 'r')
dev_st_in=open(sys.argv[3], 'r')

N=int(sys.argv[4])

def nonblank_lines(f):
	for l in f:
		line = l.rstrip()
		if line:
			yield line

exp = []
sim = []
dev_st = []

for line1 in nonblank_lines(exp_in):
    tmp = line1.split()
    exp.append(float(tmp[2]))

for line2 in nonblank_lines(sim_in):
    tmp = line2.split()
    sim.append(float(tmp[2]))

for line3 in nonblank_lines(dev_st_in):
    tmp = line3.split()
    dev_st.append(float(tmp[0]))

sum_chi=0
for i in range(0,N*N):
    if(dev_st[i]!=0):
        sum_chi+=(exp[i]-sim[i])**2/(dev_st[i]**2)


chi2=sum_chi/float(N*N)

print(chi2)
