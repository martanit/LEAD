import numpy as np
import sys
import subprocess

traj_file = sys.argv[1]
p = int(sys.argv[2])
mddt = float(sys.argv[3])
save = int(sys.argv[4])

command = "grep -c ^$ "+(traj_file)
nframe=int(subprocess.check_output(command, shell=True))

command2 = "head -n 1 "+(traj_file)
N=int(subprocess.check_output(command2, shell=True))

# Define function to read from a file without blank lines
def nonblank_lines(f):
	for l in f:
		line = l.rstrip()
		if line:
			yield line

# Define function to read xyz and divide into frame
def read_xyz(file_name):
    #Set -1 because first nonblank line of traj is always != 'Au'
    #and I need to start with the frame "0"
    frame = -1;
    m = 0;
    # Create a list containing nframe list each of N items all set to 0
    x = [[0 for t in range(int(nframe/save))] for k in range(N)]
    y = [[0 for t in range(int(nframe/save))] for k in range(N)]
    z = [[0 for t in range(int(nframe/save))] for k in range(N)]
    
    with open(file_name) as traj:
        for line in nonblank_lines(traj):
            if ( line.split()[0] == 'Au' ):
                if(frame % save == 0): 
                    name, x_c, y_c, z_c = line.split();
                    x[m][int(frame/save)] = float(x_c)
                    y[m][int(frame/save)] = float(y_c)
                    z[m][int(frame/save)] = float(z_c)
                    m = m+1
            else: 
                frame=frame+1
                m = 0
    return x,y,z

# Define function to calculate distance
def dist(x,y,z,i,j,f):
    # Call func
    return np.sqrt((x[i][f]-x[j][f])**2+(y[i][f]-y[j][f])**2+(z[i][f]-z[j][f])**2)

#################################################################################
# MAIN
x,y,z = read_xyz(traj_file)

#for i in range(int(nframe/save)):
#    print(str(N)+'\n')
#    for j in range(N):
#        print('Au'+'\t'+str(x[j][i])+'\t'+str(y[j][i])+'\t'+str(z[j][i]))



sdx = [0 for t in range(int(nframe/save))]
sdy = [0 for t in range(int(nframe/save))]
sdz = [0 for t in range(int(nframe/save))]
cnt = [0 for t in range(int(nframe/save))]

for t in range(int(nframe/save)):
    for dt in range(1,int((nframe/save)-t)):
        sdx[dt] += (x[p][t+dt] - x[p][t])**2
        sdy[dt] += (y[p][t+dt] - y[p][t])**2
        sdz[dt] += (z[p][t+dt] - z[p][t])**2
        cnt[dt] += 1
t=0
for t in range(int((nframe/save)/2)):
    if cnt[t] != 0:
        sdx[t] /= cnt[t]
        sdy[t] /= cnt[t]
        sdz[t] /= cnt[t]

    print(t*mddt, sdx[t]+sdy[t]+sdz[t])
