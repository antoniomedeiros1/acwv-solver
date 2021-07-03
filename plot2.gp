stat "data0.dat" name "A"
set xlabel "Largura (m)"
set ylabel "Profundidade (m)"
set grid
set xrange [A_min_x:A_max_x]                                              
set yrange [A_max_y:A_min_y]
set cbrange [-2:2]
load './config/magma.pal'
set output 'plot2.gif'  
set term gif animate delay 3
do for [i=0:79] {plot 'data'.i.'.dat' with image}
unset output