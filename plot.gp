cd "../mdf-onda-acustica"
stat "data0.dat" name "A"
set xlabel "x"
set ylabel "z"
set grid
set xrange [A_min_x:A_max_x]                                              
set yrange [A_min_y:A_max_y]
set cbrange [-3:3]
set palette rgb 21,22,23
set output 'plot1.gif'  
set term gif animate delay 15
do for [i=0:19] {plot 'data'.i.'.dat' with image}
unset output