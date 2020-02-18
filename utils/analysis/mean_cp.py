import numpy as np
import sys
import subprocess
lg=int(sys.argv[1])
out=sys.argv[2]
# Define function to read from a file without blank lines
def nonblank_lines(f):
	for l in f:
		line = l.rstrip()
		if line:
	 		yield line
def split(array, nrows, ncols):
    r, h = array.shape
    return (array.reshape(h//nrows, nrows, -1, ncols)
            .swapaxes(1, 2)
            .reshape(-1, nrows, ncols))

# MAIN
ps_tot = []
ps_mean = []
for i in sys.argv[3:]:
    print(i)
    hic_in = open(i, 'r')    
    ps = []
    for line in nonblank_lines(hic_in):
        tmp = line.strip()
        tmp = tmp.split(" ")
        ps.append(float(tmp[1]))
    ps_tot.append(ps)

ps_tot=np.array(ps_tot).T.tolist()
for j in range(0, lg):
    ps_mean.append(np.mean(ps_tot[j]))

out=open(out, 'w')
for t in range(0, lg):
    out.write(str(t*3200)+'\t'+str(ps_mean[t])+'\n')
out.close()
