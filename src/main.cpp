#define _USE_MATH_DEFINES
#define STENCIL 2

#include <iostream>

#include "../include/Solver2d.h"

using namespace std;

int main(int argc, char *argv[]) {
    PetscInitialize(&argc, &argv, NULL, NULL);
    if (argc == 5) {
        string input_file = argv[1];
        string output_folder = argv[2];
        int number_of_steps = atoi(argv[3]);
        float dt = atof(argv[4]);
        Solver2d solver(input_file, output_folder, number_of_steps, dt, 40);
        solver.solve();
    }
    else if (argc == 6){
        string input_file = argv[1];
        string output_folder = argv[2];
        int number_of_steps = atoi(argv[3]);
        float dt = atof(argv[4]);
        int number_of_frames = atoi(argv[5]);
        Solver2d solver(input_file, output_folder, number_of_steps, dt, number_of_frames);
        solver.solve();
    }
    else {
        string input_file = "../data/input.vti";
        string output_folder = "../data/";
        int number_of_steps = 8000;
        float dt = 0.00025;
        int number_of_frames = 10;
        Solver2d solver(input_file, output_folder, number_of_steps, dt, number_of_frames);
        solver.solve();
    }
    PetscFinalize();
    return 0;
}
