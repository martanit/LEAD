set output filename . '.svg'
set title filename
set xlabel "Log(l)"
set ylabel "Log(P(l)) []"
set terminal svg
a=1.
b=1.
f(x)=a+b*x
s(x)=c+d*x
fit [log(2000):log(100000)] f(x) filename u (log($1)):(log($2)) via a,b
fit [log(100000):log(1000000)] s(x) filename u (log($1)):(log($2)) via c,d
ti = sprintf("fit:\t%.2f", b)
tti = sprintf("fit:\t%.2f", d)
set log xy
plot [1000:1000000] filename title "simulation" w l, ((2000 < x & x<100000)?exp(f(log(x))):1/0) t ti, ((100000<x & x<1000000)?exp(s(log(x))):1/0) t tti
