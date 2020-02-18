set terminal svg
set output filename . '.svg'
set title filename
set xlabel "Log(N)"
set ylabel "Log(Rg^2) [a^2]"
a=1.
b=1.
f(x)=a+b*x
fit f(x) filename u (log($1)):(log($2)) via a,b
ti = sprintf("fit:\t %.2f", b)
set log xy
plot filename title "simulation" w l, exp(f(log(x))) t ti
