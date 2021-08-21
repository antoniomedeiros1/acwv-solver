#define _USE_MATH_DEFINES
#define STENCIL 3 

#include <iostream>
#include "omp.h"

#include "Solver3d.h"

using namespace std;

int main() {

    Solver3d solver;
    solver.solve();

    return 0;
}
