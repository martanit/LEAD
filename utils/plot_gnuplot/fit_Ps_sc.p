set fit logfile '/dev/null'
a=0.1
b=1
f(x)=a+b*x
fit [log(3000):log(100000)] f(x) filename u (log($1)):(log($2)) via a,b
ti = sprintf("%.2f", b)
set log xy
p exp(f(log(x))) t ti, filename w l
pause -1
