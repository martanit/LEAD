set fit logfile '/dev/null'
set terminal svg
set output filename . '.svg'
set title "Contact probability for Rouse polymer"
set xlabel "Log(l) [a.u.]"
set ylabel "Log(P(l)) [a.u.]"
a=1.
b=1.
f(x)=a+b*x
fit [log(2):log(350)] f(x) filename u (log($1)):(log($2)) via a,b
ti = sprintf("fit:\t%.2f", b)
set log xy
plot filename title "simulation" w l, ((2<x & x<350)?exp(f(log(x))):1/0) lw 2 t ti
