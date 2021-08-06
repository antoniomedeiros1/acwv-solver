#define _USE_MATH_DEFINES
#define STENCIL 3 

#include <iostream>
#include "Solver2d.h"

using namespace std;

int main() {

    Solver2d s("reservatorio.txt");
    s.solve();

    return 0;
}
