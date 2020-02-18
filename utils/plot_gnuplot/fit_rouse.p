set output filename . '.svg'
set title filename
set xlabel "Log(l)"
set ylabel "Log(P(l)) []"
set terminal svg
a=1.
b=1.
f(x)=a+b*x
fit [log(6000):log(1000000)] f(x) filename u (log($1)):(log($2)) via a,b
ti = sprintf("fit:\t%.2f", b)
set log xy
p exp(f(log(x))) t ti, filename title "simulation" w l
