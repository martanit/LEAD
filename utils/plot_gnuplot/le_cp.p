set fit logfile '/dev/null'
a=1.
b=1.
f(x)=a+b*x
fit [log(30000):log(300000)] f(x) filename u (log($1)):(log($2)) via a,b
ti = sprintf("%.2f", b)
set log xy
plot filename w l, ((30000<x & x<300000)?exp(f(log(x))):1/0) lw 2 t ti
pause -1
