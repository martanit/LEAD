set output filename . '.svg'
set title filename
set xlabel "Log(N)"
set ylabel "Log(Rg^2) [a^2]"
set terminal svg
a=1.
b=1.
f(x)=a+b*x
s(x)=c+d*x
fit [log(2):log(10)] f(x) filename u (log($1)):(log($2)) via a,b
fit [log(10):log(100)] s(x) filename u (log($1)):(log($2)) via c,d
ti = sprintf("fit:\t %.2f", b)
tti = sprintf("fit:\t %.2f", d)
set log xy
plot filename title "simulation" w l, ((2 < x & x<10)?exp(f(log(x))):1/0) t ti, ((10<x & x<100)?exp(s(log(x))):1/0) t tti
