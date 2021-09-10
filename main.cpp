#define _USE_MATH_DEFINES
#define STENCIL 3 

#include <iostream>
#include "omp.h"

#include "Solver3d.h"

using namespace std;

int main() {

    system("export OMP_NUM_THREADS=1");

    omp_set_num_threads(1);
	  omp_set_dynamic(0);

    int t;
    
    #pragma omp parallel
    t = omp_get_num_threads();

    cout << "\nQuantidade de threads disponiveis: " << t << "\n\n";

    //Solver3d solver;
    //solver.solve();

    Solver2d solver("double_slit.txt");
    solver.solve();

    return 0;
}
