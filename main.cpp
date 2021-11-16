#define _USE_MATH_DEFINES
#define STENCIL 3 

#include <iostream>
#include "omp.h"

#include "include/Solver3d.h"

using namespace std;

int main(int argc, char *argv[]) {

    string arq(argv[1]);

    system("export OMP_NUM_THREADS=8");

    omp_set_num_threads(8);
	omp_set_dynamic(0);

    int t;
    
    #pragma omp parallel
    t = omp_get_num_threads();

    cout << "\nQuantidade de threads disponiveis: " << t << "\n\n";

    Solver2d solver(arq);
    solver.solve();

    return 0;
}
