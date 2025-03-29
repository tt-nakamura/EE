reset
set terminal postscript eps enhanced 24
set output 'fig1.eps'
set xlabel 'length of N({/Symbol p})  / bit'
set ylabel 'time to compute \{{/Symbol a}|{/Symbol p}\}  / msec'
set key left Left reverse
plot \
	'fig1.txt' u 1:($2/0.001) t 'naive method' w l lt 2,\
	'fig1.txt' u 1:($3/0.001) t 'fast method' w l lt 1
