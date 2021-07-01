rm *.dat
cmake build
cd build
make
./mdf.exe
cd ..
gnuplot plot.gp
code plot1.gif