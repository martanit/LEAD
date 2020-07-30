#!/usr/bin/env python3

import numpy as np
import sys
import subprocess
from scipy.special import binom
from scipy.special import hyp2f1
def choose(n,k):
    if(0<= k <= n):
        return binom(n,k)
    else:
        return 0
       

ctcf_file = sys.argv[1]

command = "grep -c ^ "+(ctcf_file)
N=int(subprocess.check_output(command, shell=True))

with open(ctcf_file) as f:
    ctcf = [line.rstrip() for line in f]

P_eff = np.zeros((N,N))
U_eff = np.zeros((N,N))

kon=2.28e-15
koff=1.28e-15
kfw=2e-13
p0=kon/(koff+kfw)
for i in range(0,N-2):
    for j in range(i+2,N):
        cr=0
        cl=0
        kcr = []
        kcl = []
        for t in range(i, j+1):
            if(int(ctcf[t])>0):
                cr=cr+1
                kcr.append(t);
            if(int(ctcf[t])<0):
                cl=cl+1
                kcl.append(t);
       
        #legend: >:ic   <:jc    i<:i=jc    >j:j=ic
       
        #i!=ic and j!=jc
        ################

        #no ctcf between i and j
        #(a) i-------j
        if(cr==0 and cl==0): 
            P_eff[i,j]=(kfw/(koff+kfw))**(j-i-1)*p0
        #(b) i------>j
        if(cr==1 and cl==0):
            if(j==kcr[cr-1]):
                P_eff[i,j]=(kfw/(koff+kfw))**(j-i-1)*p0
        #(c) i<------j
        if(cl==1 and cr==0):
            if(i==kcl[0]):
                P_eff[i,j]=(kfw/(koff+kfw))**(j-i-1)*p0
        #(d) i<----->j
        if(cl==1 and cr==1):
            if(i==kcl[0] and j==kcr[0]):
                P_eff[i,j]=(kfw/(koff+kfw))**(j-i-1)*p0

        #two or more ctcf with convergent orientation
        #(e) i-->--<--j
        #(f) i-->--<->j
        #(g) i<->--<--j
        #(h) i<->--<->j
        if(cr>0 and cl>0):
           # if(i!=kcr[0] and j!=kcl[cl-1] and kcl[cl-1]!=(kcr[0]-1)):
            if(i!=kcr[0] and j!=kcl[cl-1] and kcl[cl-1]>kcr[0]):
                P_eff[i,j]=0

        #one or several same oriented ctcf between i and j
        if(cr>0 and cl==0):
            #(i) i-->-->--j
            if(i!=kcr[0] and j!=kcr[cr-1]):
                n=kcr[0]-i
                m=j-kcr[cr-1]
                P_eff[i,j]=choose(n+m-1,m)*hyp2f1(1,1-n,1+m,-1)*(kfw/(2*(koff+kfw)))**(j-i-1)*p0
        if(cr>1 and cl==0):        
            #(k) i-->-->->j
            if(i!=kcr[0] and j==kcr[cr-1]):
                n=kcr[0]-i
                m=j-kcr[cr-2]
                P_eff[i,j]=choose(n+m-1,m)*hyp2f1(1,1-n,1+m,-1)*(kfw/(2*(koff+kfw)))**(j-i-1)*p0
        if(cl>0 and cr==0):
            #(j) i--<--<--j
            if(i!=kcl[0] and j!=kcl[cl-1]):
                m=kcl[0]-i
                n=j-kcl[cl-1]
                P_eff[i,j]=choose(n+m-1,m)*hyp2f1(1,1-n,1+m,-1)*(kfw/(2*(koff+kfw)))**(j-i-1)*p0
        if(cl>1 and cr==0):        
            #(l) i<--<--<--j
            if(i==kcl[0] and j!=kcl[cl-1]):
                m=kcl[1]-i
                n=j-kcl[cl-1]
                P_eff[i,j]=choose(n+m-1,m)*hyp2f1(1,1-n,1+m,-1)*(kfw/(2*(koff+kfw)))**(j-i-1)*p0
        if(cr>0 and cl==1):
            #(m) i<->-->--j
            if(i==kcl[0] and j!=kcr[cr-1]):
                n=kcr[0]-i
                m=j-kcr[cr-1]
                P_eff[i,j]=choose(n+m-1,m)*hyp2f1(1,1-n,1+m,-1)*(kfw/(2*(koff+kfw)))**(j-i-1)*p0        
        if(cl>1 and cr==1):    
            #(p) i<-->-->-->j
            if(i==kcl[0] and j==kcr[cr-1]):
                n=kcr[0]-i
                m=j-kcl[cl-2]
                P_eff[i,j]=choose(n+m-1,m)*hyp2f1(1,1-n,1+m,-1)*(kfw/(2*(koff+kfw)))**(j-i-1)*p0
        if(cl>0 and cr==1):
            #(n) i--<--<->j
            if(j==kcr[cr-1] and i!=kcl[0]):
                m=kcl[0]-i
                n=j-kcl[cl-1]
                P_eff[i,j]=choose(n+m-1,m)*hyp2f1(1,1-n,1+m,-1)*(kfw/(2*(koff+kfw)))**(j-i-1)*p0
        if(cl>1 and cr==1):        
            #(o) i<--<--<-->j
            if(i==kcl[0] and j==kcr[cr-1]):
                m=kcl[1]-i
                n=j-kcl[cl-1]
                P_eff[i,j]=choose(n+m-1,m)*hyp2f1(1,1-n,1+m,-1)*(kfw/(2*(koff+kfw)))**(j-i-1)*p0
        
        #two or more ctcf with divergent orientation
        if(cr>0 and cl>0):
            #(q) i---<<-->>---j
            if(i!=kcr[0] and j!=kcl[cl-1] and kcl[cl-1]<kcr[0] and i!=kcl[0] and j!=kcr[cr-1]):
                n=kcl[0]-i
                m=j-kcr[cr-1]
                P_eff[i,j]=((kfw/2.)/(koff+kfw))**(j-i-1)*(choose(m-1,m)*hyp2f1(1,1,1+m,-1)*((kfw/2.)/(koff+kfw/2.))**n+\
                                                           choose(n-1,0)*hyp2f1(1,1-n,1,-1)*((kfw/2.)/(koff+kfw/2.))**m)*p0
        if(cr>0 and cl>1):
            #(r) i<--<<-->>---j   
            if(i!=kcr[0] and j!=kcl[cl-1] and kcl[cl-1]<kcr[0] and i==kcl[0] and j!=kcr[cr-1]):
                n=kcl[1]-i
                m=j-kcr[cr-1]
                P_eff[i,j]=((kfw/2.)/(koff+kfw))**(j-i-1)*(choose(m-1,m)*hyp2f1(1,1,1+m,-1)*((kfw/2.)/(koff+kfw/2.))**n+\
                                                           choose(n-1,0)*hyp2f1(1,1-n,1,-1)*((kfw/2.)/(koff+kfw/2.))**m)*p0
        if(cr>1 and cl>0):
            #(s) i---<<-->>-->j
            if(i!=kcr[0] and j!=kcl[cl-1] and kcl[cl-1]<kcr[0] and i!=kcl[0] and j==kcr[cr-1]):
                n=kcl[0]-i
                m=j-kcr[cr-2]
                P_eff[i,j]=((kfw/2.)/(koff+kfw))**(j-i-1)*(choose(m-1,m)*hyp2f1(1,1,1+m,-1)*((kfw/2.)/(koff+kfw/2.))**n+\
                                                           choose(n-1,0)*hyp2f1(1,1-n,1,-1)*((kfw/2.)/(koff+kfw/2.))**m)*p0
        if(cr>1 and cl>1):
            #(t) i<--<<-->>-->j
            if(i!=kcr[0] and j!=kcl[cl-1] and kcl[cl-1]<kcr[0] and i==kcl[0] and j==kcr[cr-1]):
                n=kcl[1]-i
                m=j-kcr[cr-2]
                P_eff[i,j]=((kfw/2.)/(koff+kfw))**(j-i-1)*(choose(m-1,m)*hyp2f1(1,1,1+m,-1)*((kfw/2.)/(koff+kfw/2.))**n+\
                                                           choose(n-1,0)*hyp2f1(1,1-n,1,-1)*((kfw/2.)/(koff+kfw/2.))**m)*p0
        
        #i=ic and j!=jc
        ################
        
        #no ctcf between i and j           
        if((cr==1 or cr==2) and cl==0 ):
            #(u) i>------j
            if(cr==1 and i==kcr[0]):
                P_eff[i,j]=kon/(koff+kfw/2.)*((kfw/2.)/(koff+kfw/2.))**(j-i-2)\
                            +(koff+kfw)/(kfw/2.)*((kfw/2.)/(koff+kfw/2.))**(j-i-2)*p0-((kfw/2.)/(koff+kfw))**(j-i-3)*p0
            #(v) i>----->j
            if(cr==2 and i==kcr[0] and j==kcr[cr-1]):
                P_eff[i,j]=kon/(koff+kfw/2.)*((kfw/2.)/(koff+kfw/2.))**(j-i-2)\
                            +(koff+kfw)/(kfw/2.)*((kfw/2.)/(koff+kfw/2.))**(j-i-2)*p0-((kfw/2.)/(koff+kfw))**(j-i-3)*p0
        
        #two or more ctcf with convergent orientation
        #(w) i>-->--<--j
        #(x) i>-->--<->j
        if(cr>1 and cl>0):
            if(i==kcr[0] and j!=kcl[cl-1] and kcl[cl-1]>kcr[1]):
                P_eff[i,j]=0

        #one or several same oriented ctcf between i and j
        if(cr>1 and cl==0):
            #(y) i>-->-->--j
            if(i==kcr[0] and j!=kcr[cr-1]):
                n1=kcr[1]-i-1
                m1=j-kcr[cr-1]
                n2=kcr[1]-i
                m2=j-kcr[cr-1]-1
                P_eff[i,j]=(kfw/2.)/(koff+kfw/2.)*(choose(n1+m1-1,m1)*hyp2f1(1,1-n1,1+m1,-1)+\
                                                   choose(n2+m2-1,m2)*hyp2f1(1,1-n2,1+m2,-1))*(kfw/(2*(koff+kfw)))**(j-i-2)*p0
        if(cr>2 and cl==0):
            #(a') i>-->-->->j
            if(i==kcr[0] and j==kcr[cr-1]):
                n1=kcr[1]-i-1
                m1=j-kcr[cr-2]
                n2=kcr[1]-i
                m2=j-kcr[cr-2]-1
                P_eff[i,j]=(kfw/2.)/(koff+kfw/2.)*(choose(n1+m1-1,m1)*hyp2f1(1,1-n1,1+m1,-1)+\
                                                   choose(n2+m2-1,m2)*hyp2f1(1,1-n2,1+m2,-1))*(kfw/(2*(koff+kfw)))**(j-i-2)*p0
        if(cl>0 and cr==1):
            #(z) i>--<--<--j
            if(i==kcr[0] and j!=kcl[cl-1] and j!=kcr[cr-1]):
                m1=kcl[0]-i-1
                n1=j-kcl[cl-1]
                m2=kcl[0]-i
                n2=j-kcl[cl-1]-1
                P_eff[i,j]=(kfw/2.)/(koff+kfw/2.)*(choose(n1+m1-1,m1)*hyp2f1(1,1-n1,1+m1,-1)+\
                                                   choose(n2+m2-1,m2)*hyp2f1(1,1-n2,1+m2,-1))*(kfw/(2*(koff+kfw)))**(j-i-2)*p0
        if(cl>0 and cr==2):
            #(b') i>--<--<->j         
            if(i==kcr[0] and j==kcr[cr-1]):
                m1=kcl[0]-i-1
                n1=j-kcl[cl-1]
                m2=kcl[0]-i
                n2=j-kcl[cl-1]-1
                P_eff[i,j]=(kfw/2.)/(koff+kfw/2.)*(choose(n1+m1-1,m1)*hyp2f1(1,1-n1,1+m1,-1)+\
                                                   choose(n2+m2-1,m2)*hyp2f1(1,1-n2,1+m2,-1))*(kfw/(2*(koff+kfw)))**(j-i-2)*p0
        
        #two or more ctcf with divergent orientation
        if(cr>1 and cl>0):
            #(c') i>--<<-->>--j
            if(i==kcr[0] and j!=kcl[cl-1] and j!=kcr[cr-1] and kcl[cl-1]<kcr[1]):
                n1=kcl[0]-i-1
                m1=j-kcr[cr-1]
                n2=kcl[0]-i
                m2=j-kcr[cr-1]-1
                P_eff[i,j]=(kfw/2.)/(koff+kfw/2.)*(choose(m1-1,m1)*hyp2f1(1,1,1+m1,-1)*((kfw/2.)/(koff+kfw))**(m1+kcr[cr-1]-kcl[0]-1)*((kfw/2.)/(koff+kfw/2.))**n1+\
                                          choose(n1-1,0)*hyp2f1(1,1-n1,1,-1)*((kfw/2.)/(koff+kfw))**(n1+kcr[cr-1]-kcl[0]-1)*((kfw/2.)/(koff+kfw/2.))**m1+\
                                          choose(m2-1,m2)*hyp2f1(1,1,1+m2,-1)*((kfw/2.)/(koff+kfw))**(m2+kcr[cr-1]-kcl[0]-1)*((kfw/2.)/(koff+kfw/2.))**n2+\
                                          choose(n2-1,0)*hyp2f1(1,1-n2,1,-1)*((kfw/2.)/(koff+kfw))**(n2+kcr[cr-1]-kcl[0]-1)*((kfw/2.)/(koff+kfw/2.))**m2)*p0
        if(cr>2 and cl>0):
            #(d') i>--<<-->>->j
            if(i==kcr[0] and j==kcr[cr-1] and kcl[cl-1]<kcr[1]):
                n1=kcl[0]-i-1
                m1=j-kcr[cr-2]
                n2=kcl[0]-i
                m2=j-kcr[cr-2]-1
                P_eff[i,j]=(kfw/2.)/(koff+kfw/2.)*(choose(m1-1,m1)*hyp2f1(1,1,1+m1,-1)*((kfw/2.)/(koff+kfw))**(m1+kcr[cr-1]-kcl[0]-1)*((kfw/2.)/(koff+kfw/2.))**n1+\
                                          choose(n1-1,0)*hyp2f1(1,1-n1,1,-1)*((kfw/2.)/(koff+kfw))**(n1+kcr[cr-1]-kcl[0]-1)*((kfw/2.)/(koff+kfw/2.))**m1+\
                                          choose(m2-1,m2)*hyp2f1(1,1,1+m2,-1)*((kfw/2.)/(koff+kfw))**(m2+kcr[cr-1]-kcl[0]-1)*((kfw/2.)/(koff+kfw/2.))**n2+\
                                          choose(n2-1,0)*hyp2f1(1,1-n2,1,-1)*((kfw/2.)/(koff+kfw))**(n2+kcr[cr-1]-kcl[0]-1)*((kfw/2.)/(koff+kfw/2.))**m2)*p0

        #i!=ic and j=jc
        ################

        #no ctcf between i and j           
        if(cr==0 and (cl==1 or cl==2)):
            #(e') i------<j
            if(cl==1 and j==kcl[cl-1]):
                P_eff[i,j]=kon/(koff+kfw/2.)*((kfw/2.)/(koff+kfw/2.))**(j-i-2)\
                                +(koff+kfw)/(kfw/2.)*((kfw/2.)/(koff+kfw/2.))**(j-i-2)*p0-((kfw/2.)/(koff+kfw))**(j-i-3)*p0
            #(f') i<------<j
            if(cl==2 and i==kcl[0] and j==kcl[cl-1]):
                P_eff[i,j]=kon/(koff+kfw/2.)*((kfw/2.)/(koff+kfw/2.))**(j-i-2)\
                                +(koff+kfw)/(kfw/2.)*((kfw/2.)/(koff+kfw/2.))**(j-i-2)*p0-((kfw/2.)/(koff+kfw))**(j-i-3)*p0

        #two or more ctcf with convergent orientation
        #(g') i-->--<--<j
        #(h') i<->--<--<j
        if(cr>0 and cl>1):
            if(i!=kcr[0] and j==kcl[cl-1] and kcl[cl-2]!=(kcr[0]-1)):
                P_eff[i,j]=0

        #one or several same oriented ctcf between i and j
        if(cl>1 and cr==0):
            #(i') i--<--<-<j
            if(j==kcl[cl-1] and i!=kcl[0]):
                m1=kcl[0]-i-1
                n1=j-kcl[cl-2]
                m2=kcl[0]-i
                n2=j-kcl[cl-2]-1
                P_eff[i,j]=(kfw/2.)/(koff+kfw/2.)*(choose(n1+m1-1,m1)*hyp2f1(1,1-n1,1+m1,-1)+\
                                                   choose(n2+m2-1,m2)*hyp2f1(1,1-n2,1+m2,-1))*(kfw/(2*(koff+kfw)))**(j-i-2)*p0
        if(cl>2 and cr==0):
            #(k') i<-<--<-<j
            if(j==kcl[cl-1] and i==kcl[0]):
                m1=kcl[1]-i-1
                n1=j-kcl[cl-2]
                m2=kcl[1]-i
                n2=j-kcl[cl-2]-1
                P_eff[i,j]=(kfw/2.)/(koff+kfw/2.)*(choose(n1+m1-1,m1)*hyp2f1(1,1-n1,1+m1,-1)+\
                                                   choose(n2+m2-1,m2)*hyp2f1(1,1-n2,1+m2,-1))*(kfw/(2*(koff+kfw)))**(j-i-2)*p0
        #(j') i-->-->-<j
        if(cr>0 and cl==1):
            if(j==kcl[cl-1] and i!=kcr[0]):
                n1=kcr[0]-i-1
                m1=j-kcr[cr-1]
                n2=kcr[0]-i
                m2=j-kcr[cr-1]-1
                P_eff[i,j]=(kfw/2.)/(koff+kfw/2.)*(choose(n1+m1-1,m1)*hyp2f1(1,1-n1,1+m1,-1)+\
                                                   choose(n2+m2-1,m2)*hyp2f1(1,1-n2,1+m2,-1))*(kfw/(2*(koff+kfw)))**(j-i-2)*p0
        #(l') i<->-->-<j       
        if(cr>0 and cl==2):
            if(i==kcl[0] and j==kcl[cl-1]):
                n1=kcr[0]-i-1
                m1=j-kcr[cr-1]
                n2=kcr[0]-i
                m2=j-kcr[cr-1]-1
                P_eff[i,j]=(kfw/2.)/(koff+kfw/2.)*(choose(n1+m1-1,m1)*hyp2f1(1,1-n1,1+m1,-1)+\
                                                   choose(n2+m2-1,m2)*hyp2f1(1,1-n2,1+m2,-1))*(kfw/(2*(koff+kfw)))**(j-i-2)*p0
        
        #two or more ctcf with divergent orientation
        if(cr>0 and cl>1):
            #(m') i--<<-->>-<j
            if(i!=kcr[0] and j==kcl[cl-1] and i!=kcl[0] and kcl[cl-2]<kcr[0]):
                n1=kcl[0]-i-1
                m1=j-kcr[cr-1]
                n2=kcl[0]-i
                m2=j-kcr[cr-1]-1
                P_eff[i,j]=(kfw/2.)/(koff+kfw/2.)*(choose(m1-1,m1)*hyp2f1(1,1,1+m1,-1)*((kfw/2.)/(koff+kfw))**(m1+kcr[cr-1]-kcl[0]-1)*((kfw/2.)/(koff+kfw/2.))**n1+\
                                          choose(n1-1,0)*hyp2f1(1,1-n1,1,-1)*((kfw/2.)/(koff+kfw))**(n1+kcr[cr-1]-kcl[0]-1)*((kfw/2.)/(koff+kfw/2.))**m1+\
                                          choose(m2-1,m2)*hyp2f1(1,1,1+m2,-1)*((kfw/2.)/(koff+kfw))**(m2+kcr[cr-1]-kcl[0]-1)*((kfw/2.)/(koff+kfw/2.))**n2+\
                                          choose(n2-1,0)*hyp2f1(1,1-n2,1,-1)*((kfw/2.)/(koff+kfw))**(n2+kcr[cr-1]-kcl[0]-1)*((kfw/2.)/(koff+kfw/2.))**m2)*p0
        if(cr>0 and cl>2):
            #(n') i<--<<-->>-<j
            if(i==kcl[0] and j==kcl[cl-1] and kcl[cl-2]<kcr[0]):
                n1=kcl[1]-i-1
                m1=j-kcr[cr-1]
                n2=kcl[1]-i
                m2=j-kcr[cr-1]-1
                P_eff[i,j]=(kfw/2.)/(koff+kfw/2.)*(choose(m1-1,m1)*hyp2f1(1,1,1+m1,-1)*((kfw/2.)/(koff+kfw))**(m1+kcr[cr-1]-kcl[0]-1)*((kfw/2.)/(koff+kfw/2.))**n1+\
                                          choose(n1-1,0)*hyp2f1(1,1-n1,1,-1)*((kfw/2.)/(koff+kfw))**(n1+kcr[cr-1]-kcl[0]-1)*((kfw/2.)/(koff+kfw/2.))**m1+\
                                          choose(m2-1,m2)*hyp2f1(1,1,1+m2,-1)*((kfw/2.)/(koff+kfw))**(m2+kcr[cr-1]-kcl[0]-1)*((kfw/2.)/(koff+kfw/2.))**n2+\
                                          choose(n2-1,0)*hyp2f1(1,1-n2,1,-1)*((kfw/2.)/(koff+kfw))**(n2+kcr[cr-1]-kcl[0]-1)*((kfw/2.)/(koff+kfw/2.))**m2)*p0
        
        #i=ic and j=jc
        ################
        
        #no ctcf between i and j
        #(o') i>------<j
        if(cr==1 and cl==1):
            if(i==kcr[0] and j==kcl[cl-1]):
                P_eff[i,j]=(kfw/koff)*(kon/(koff+kfw/2.)*((kfw/2.)/(koff+kfw/2.))**(j-i-2)\
                                    +(koff+kfw)/(kfw/2.)*((kfw/2.)/(koff+kfw/2.))**(j-i-2)*p0-((kfw/2.)/(koff+kfw))**(j-i-3)*p0)
        #two or more ctcf with convergent orientation
        #(p') i>->--<-<j
        if(cr>1 and cl>1):
            if(i==kcr[0] and j==kcl[cl-1] and kcl[cl-2]>kcr[1]):
                P_eff[i,j]=0
        
        #one or several same oriented ctcf between i and j
        #(q') i>->-->-<j
        if(cr>1 and cl==1):
            if(i==kcr[0] and j==kcl[cl-1]):
                n1=kcr[1]-i-1
                m1=j-kcr[cr-1]
                n2=kcr[1]-i
                m2=j-kcr[cr-1]-1
                P_eff[i,j]=(kfw/2.)/koff*(choose(n1+m1-1,m1)*hyp2f1(1,1-n1,1+m1,-1)+choose(n2+m2-1,m2)*hyp2f1(1,1-n2,1+m2,-1))*(kfw/(2*(koff+kfw)))**(j-i-2)*p0
        
        #(r') i>-<--<-<j
        if(cl>1 and cr==1):
            if(i==kcr[0] and j==kcl[cl-1]):
                m1=kcl[0]-i-1
                n1=j-kcl[cl-2]
                m2=kcl[0]-i
                n2=j-kcl[cl-2]-1
                P_eff[i,j]=(kfw/2.)/koff*(choose(n1+m1-1,m1)*hyp2f1(1,1-n1,1+m1,-1)+choose(n2+m2-1,m2)*hyp2f1(1,1-n2,1+m2,-1))*(kfw/(2*(koff+kfw)))**(j-i-2)*p0
        
        #two or more ctcf with divergent orientation
        if(cr>1 and cl>1):
            #(s') i>--<<-->>--<j
            if(i==kcr[0] and j==kcl[cl-1] and kcl[cl-2]<kcr[1]):
                n1=kcl[0]-i-1
                m1=j-kcr[cr-1]
                n2=kcl[0]-i
                m2=j-kcr[cr-1]-1
                P_eff[i,j]=(kfw/2.)/koff*(choose(m1-1,m1)*hyp2f1(1,1,1+m1,-1)*((kfw/2.)/(koff+kfw))**(m1+kcr[cr-1]-kcl[0]-1)*((kfw/2.)/(koff+kfw/2.))**n1+\
                                          choose(n1-1,0)*hyp2f1(1,1-n1,1,-1)*((kfw/2.)/(koff+kfw))**(n1+kcr[cr-1]-kcl[0]-1)*((kfw/2.)/(koff+kfw/2.))**m1+\
                                          choose(m2-1,m2)*hyp2f1(1,1,1+m2,-1)*((kfw/2.)/(koff+kfw))**(m2+kcr[cr-1]-kcl[0]-1)*((kfw/2.)/(koff+kfw/2.))**n2+\
                                          choose(n2-1,0)*hyp2f1(1,1-n2,1,-1)*((kfw/2.)/(koff+kfw))**(n2+kcr[cr-1]-kcl[0]-1)*((kfw/2.)/(koff+kfw/2.))**m2)*p0
                
                
with open("prob_eff.dat",'w') as pr_out:
    for i in range(0,N-2):
        pr_out.write('\n')
        for j in range(i+2,N):
            pr_out.write(str(i)+" "+str(j)+" "+str(P_eff[i,j])+'\n')

Umax=-9999
for i in range(0,N-2):
    for j in range(i+2,N):
       if P_eff[i,j]>0:
            U_eff[i,j]=-np.log(P_eff[i,j])
       else:
            U_eff[i,j]=-1
       if U_eff[i,j]>Umax: Umax=U_eff[i,j]

#U_eff_max=U_eff[np.isfinite(U_eff)].max()

with open("pot_eff.in",'w') as po_out:
    for i in range(0,N-2):
        for j in range(i+2,N):
           if U_eff[i,j]>0:
                po_out.write(str(i)+" "+str(j)+" "+str(U_eff[i][j]-Umax)+'\n')
           else:
                po_out.write(str(i)+" "+str(j)+" "+str(0)+'\n')

