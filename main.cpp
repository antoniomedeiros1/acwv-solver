#define _USE_MATH_DEFINES
#define STENCIL 3 

#include <iostream>
#include "omp.h"

#include "Solver3d.h"

using namespace std;

int main() {

    Solver2d solver("semicirculo.txt");
    solver.solve();

    return 0;
}
