a=0.0076
b=0.9999
c=1
f(x)=a*c*b**x
fit [0:15000] f(x) filename  via b,c
plot [0:15000] filename w l, f(x)
pause -1
