set fit logfile '/dev/null'
set output filename . '.svg'
set title filename
set xlabel "tau [s]"
set ylabel "msd [a^2]"
set terminal svg
a=1
b=1
f(x)=a+b*sqrt(x)
s(x)=c+d*x
fit [0:100] f(x) filename u 1:2 via a, b
fit [100:500] s(x) filename u 1:2 via c, d
ti = sprintf("fit:\t%.2f+%.2fsqrt(x)", a, b)
si = sprintf("fit:\t%.2f+%.2fx", c, d)
p [0:500] f(x) t ti, s(x) t si, filename u 1:2 title "simulation" w l
