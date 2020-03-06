import sys
import numpy as np

repl1=open(sys.argv[1], 'r')
repl2=open(sys.argv[2], 'r')

# Define function to read from a file without blank lines
def nonblank_lines(f):
	for l in f:
		line = l.rstrip()
		if line:
			yield line

contact1 = []
contact2 = []

for line1 in nonblank_lines(repl1):
    tmp = line1.split()
    contact1.append(float(tmp[2]))

for line2 in nonblank_lines(repl2):
    tmp = line2.split()
    contact2.append(float(tmp[2]))

a = np.array([contact1, contact2])
std=np.std(a, axis=0)

for i in range(0, len(std)):
    print(std[i])
