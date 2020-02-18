from scipy.stats import chisquare
import numpy as np
import sys
import random

# Define function to read from a file without blank lines
def nonblank_lines(f):
	for l in f:
		line = l.rstrip()
		if line:
	 		yield line
lg=103
hic_sim=[]
hic_lj = []
hic_exp=[]
hic_rnd = [random.uniform(0,1) for i in range(lg**2)]

with open(sys.argv[1]) as sim_le:
    for line in nonblank_lines(sim_le):
        hic_sim.append(float(line.split()[2]))

with open(sys.argv[2]) as sim_lj:
    for line in nonblank_lines(sim_lj):
        hic_lj.append(float(line.split()[2]))

with open(sys.argv[3]) as exp_in:
    for line in nonblank_lines(exp_in):
        hic_exp.append(float(line.split()[2]))

sum_s=0
sum_r=0
sum_lj=0

for k in range(lg):
    sum_s+=((hic_sim[k]-hic_exp[k])**2)
    sum_lj+=((hic_lj[k]-hic_exp[k])**2)
    sum_r+=((hic_rnd[k]-hic_exp[k])**2)

#print(sum_s, sum_lj, sum_r)

from scipy.stats import entropy

print(entropy(hic_sim, qk=hic_exp))

