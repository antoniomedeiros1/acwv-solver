#define _USE_MATH_DEFINES
#define STENCIL 3 

#include <iostream>
#include <mpi.h>

#include "include/Solver2d.h"

using namespace std;

int main(int argc, char *argv[]) {
    if (argc == 5) {
        string input_file = argv[1];
        string output_file = argv[2];
        int number_of_steps = atoi(argv[3]);
        int dt = atoi(argv[4]);
        printf("Input file: %s\n", input_file.c_str());
        printf("Output file: %s\n", output_file.c_str());
        printf("Number of steps: %d\n", number_of_steps);
        printf("dt: %d\n", dt);
        Solver2d solver(input_file, output_file, number_of_steps, dt, 40);
        solver.solve();
        return 0;
    }
    else if (argc == 6){
        string input_file = argv[1];
        string output_file = argv[2];
        int number_of_steps = atoi(argv[3]);
        int dt = atoi(argv[4]);
        int number_of_frames = atoi(argv[5]);
        printf("Input file: %s\n", input_file.c_str());
        printf("Output file: %s\n", output_file.c_str());
        printf("Number of steps: %d\n", number_of_steps);
        printf("dt: %d\n", dt);
        Solver2d solver(input_file, output_file, number_of_steps, dt, number_of_frames);
        solver.solve();
        return 0;
    }
    else {
        printf("Usage: %s <input_file> <output_file> <number_of_steps> <dt> [number_of_frames]\n", argv[0]);
        return 1;
    }
}
