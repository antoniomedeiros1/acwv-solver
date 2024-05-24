rm *.dat
rm *.vti
cmake build
cd build
make
./acwv ../data/input.vti ../data/ 8000 0.00025 
cd ..