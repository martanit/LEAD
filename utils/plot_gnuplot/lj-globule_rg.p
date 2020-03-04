set output filename . '.svg'
set terminal svg
set title "Radius of gyration for Lennard Jones polymer"
set xlabel "Log(N) [a.u.]"
set ylabel "Log(Rg^2) [a^2]"
a=1.
b=1.
f(x)=a+b*x
s(x)=c+d*x
fit [log(2):log(20)] f(x) filename u (log($1)):(log($2)) via a,b
fit [log(20):log(100)] s(x) filename u (log($1)):(log($2)) via c,d
ti = sprintf("fit:\t %.2f", b)
tti = sprintf("fit:\t %.2f", d)
set key left top
set log xy
plot filename title "simulation" w l, ((2 < x & x<20)?exp(f(log(x))):1/0) t ti, ((20<x & x<100)?exp(s(log(x))):1/0) t tti
