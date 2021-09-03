#define _USE_MATH_DEFINES
#define STENCIL 3 

#include <iostream>
#include "omp.h"

#include "Solver3d.h"

using namespace std;

int main() {

    

    system("export OMP_NUM_THREADS=8");

    omp_set_num_threads(8);
	omp_set_dynamic(0);
	omp_set_schedule(omp_sched_dynamic, 1);

    int t;
    
    #pragma omp parallel
    t = omp_get_num_threads();

    cout << "\nQuantidade de threads disponiveis: " << t << "\n\n";

    // Solver3d solver;
    // solver.solve();

    Solver2d solver("semicirculo.txt");
    solver.solve();

    return 0;
}
