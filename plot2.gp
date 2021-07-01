stat "data0.dat" name "A"
set xlabel "Largura (m)"
set ylabel "Profundidade (m)"
set grid
set xrange [A_min_x:A_max_x]                                              
set yrange [A_max_y:A_min_y]
set cbrange [-1:1]
load './config/spectral.pal'
set output 'plot2.gif'  
set term gif animate delay 10
do for [i=0:39] {plot 'data'.i.'.dat' with image}
unset output