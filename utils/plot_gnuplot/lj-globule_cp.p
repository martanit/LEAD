set fit logfile '/dev/null'
set output filename . '.svg'
set title "Contact probability for Lennard Jones Globule"
set xlabel "Log(l) [a.u.]"
set ylabel "Log(P(l)) [a.u.]"
set terminal svg
a=1.
b=1.
f(x)=a+b*x
s(x)=c+d*x
fit [log(1):log(200)] f(x) filename u (log($1)):(log($2)) via a,b
fit [log(200):log(1000)] s(x) filename u (log($1)):(log($2)) via c,d
ti = sprintf("fit:\t%.2f", b)
tti = sprintf("fit:\t%.2f", d)
set log xy
plot [1:1000] filename title "simulation" w l, ((1 < x & x<200)?exp(f(log(x))):1/0) t ti, ((200<x & x<1000)?exp(s(log(x))):1/0) t tti
