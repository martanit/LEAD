import sys
import numpy as np

lg=int(sys.argv[2])
fin=open(sys.argv[1], 'r')


matrix = np.zeros((lg,lg))
for line in fin:
    tmp = line.split()
    idx = int(tmp[1])
    idy = int(tmp[2])
    if((idx <lg and idx >=0) and (idy < lg and idy >= 0)):
        matrix[idx, idy] = tmp[3]

max_el=np.amax(matrix)
out=open("norm_exp_"+str(lg)+"_hi-c_repl1.dat", 'w')
for t in range(lg):
    count=0
    for s in range(lg):
        if(count==0):
            out.write('\n')
        out.write(str(t)+'\t'+str(s)+'\t'+str(matrix[t][s]/max_el)+'\n')
        count=+1
out.close()
