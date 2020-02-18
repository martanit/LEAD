import numpy as np
import sys
import subprocess
lg=int(sys.argv[1])
cg=int(sys.argv[2])
out=sys.argv[3]
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
index_i = []
index_j = []
cg_tot = []
matrix = np.zeros((lg,lg))
for i in sys.argv[4:]:
    print(i)
    hic_in = open(i, 'r')    
    cg_mean = []
    for line in nonblank_lines(hic_in):
        tmp = line.strip()
        tmp = tmp.split(" ")
        idx = tmp[0]
        idy = tmp[1]
        matrix[int(idx), int(idy)] = tmp[2]
    for k in range(0, int(lg/cg)**2):
        cg_mean.append(np.mean(split(matrix,cg,cg)[k]))
    cg_tot.append(cg_mean)
s_mean = np.mean(cg_tot, axis=0)

h = []
for l in range(0, int(lg/cg)**2):
    a = np.zeros((cg,cg))
    for m in range (0, cg):
        for k in range(0, cg):
            a[m][k] = s_mean[l]
    h.append(a)
h=np.array(h)
h=h.reshape(int(lg/cg),int(lg/cg),cg,cg)
h=h.transpose(0,2,1,3).reshape(-1,lg)
print(h)

out=open(out, 'w')
for t in range(0, lg):
    count=0
    for s in range(0, lg):
        if(count==0):
            out.write('\n')
        out.write(str(t)+'\t'+str(s)+'\t'+str(h[t][s])+'\n')
        count=+1
out.close()
