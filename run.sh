cmake build
cd build
make
./mdf.exe
cd ..
gnuplot plot2.gp
rm *.dat
code plot2.gif