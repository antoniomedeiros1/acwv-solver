#define _USE_MATH_DEFINES
#define STENCIL 3 

#include <iostream>

#ifdef _OPENMP
#include "omp.h"
#endif

#include "include/Solver3d.h"

using namespace std;

int main(int argc, char *argv[]) {

    string arq(argv[1]);

    //system("export OMP_NUM_THREADS=8");

    int t=1;
    //omp_set_num_threads(8);
#ifdef _OPENMP
    #pragma omp parallel
    t = omp_get_num_threads();
#endif

    cout << "\nQuantidade de threads disponiveis: " << t << "\n\n";

    Solver2d solver(arq);
    solver.solve();

    return 0;
}
