set fit logfile '/dev/null'
set terminal svg
set output filename . '.svg'
set title "Mean squared displacement for Rouse polymer"
set xlabel "tau [s]"
set ylabel "msd [a^2]"
a=1
b=1
f(x)=a+b*sqrt(x)
s(x)=c+d*x
fit [0:40] f(x) filename u 1:2 via a, b
fit [40:5000] s(x) filename u 1:2 via c, d
ti = sprintf("fit:\t%.2f+%.2fsqrt(x)", a, b)
si = sprintf("fit:\t%.2f+%.2fx", c, d)
p [0:5000] f(x) t ti, s(x) t si, filename u 1:2 title "simulation" w l
