#define _USE_MATH_DEFINES
#define STENCIL 3 

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
        // printf("Input file: %s\n", input_file.c_str());
        // printf("Output file: %s\n", output_folder.c_str());
        // printf("Number of steps: %d\n", number_of_steps);
        // printf("dt: %f\n", dt);
        Solver2d solver(input_file, output_folder, number_of_steps, dt, 40);
        solver.solve();
    }
    else if (argc == 6){
        string input_file = argv[1];
        string output_folder = argv[2];
        int number_of_steps = atoi(argv[3]);
        float dt = atof(argv[4]);
        int number_of_frames = atoi(argv[5]);
        // printf("Input file: %s\n", input_file.c_str());
        // printf("Output file: %s\n", output_folder.c_str());
        // printf("Number of steps: %d\n", number_of_steps);
        // printf("dt: %f\n", dt);
        Solver2d solver(input_file, output_folder, number_of_steps, dt, number_of_frames);
        // solver.solve();
    }
    else {
        string input_file = "../data/input.vti";
        string output_folder = "../data/";
        int number_of_steps = 8000;
        float dt = 0.00025;
        int number_of_frames = 40;
        Solver2d solver(input_file, output_folder, number_of_steps, dt, number_of_frames);
        solver.solve();
        // printf("Usage: %s <input_file> <output_file> <number_of_steps> <dt> [number_of_frames]\n", argv[0]);
    }
    PetscFinalize();
    return 0;
}
