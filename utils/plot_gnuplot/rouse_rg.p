set terminal svg
set output filename . '.svg'
set title "Radius of gyration for Rouse polymer"
set xlabel "Log(N)"
set ylabel "Log(Rg^2) [a^2]"
a=1.
b=1.
f(x)=a+b*x
fit [log(2):log(100)] f(x) filename u (log($1)):(log($2)) via a,b
ti = sprintf("fit:\t %.2f", b)
set log xy
plot filename title "simulation" w l , ((2<x & x<100)?exp(f(log(x))):1/0) lw 2 title ti
