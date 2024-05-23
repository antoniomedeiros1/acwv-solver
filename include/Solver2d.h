#include <iostream>
#include <fstream>
#include <cmath>
#include "math.h"
#include <string>
#include <chrono>

#include "Domain.h"

#define STENCIL 3

using namespace std;

class Solver2d{

    public:
        Solver2d(string input_file, string output_file, int number_of_steps, int dt, int number_of_frames);
        ~Solver2d();
        void printParameters();
        void saveVTI(Domain d, Grid2d* u, string nomeDoArq, string info);
        void saveVTIbin(Domain d, Grid2d* u, string nomeDoArq, string info);
        void solve();
    
    private: 
        void readInput(string input_file);
        void computeNext(Domain d, Grid2d* u_current, Grid2d* u_next, int k);
        void computeNext(Domain d, Grid2d* u_current, Grid2d* u_next);
        float source(int x, int z, float k);
        void applyReynoldsBC(Grid2d* u_current, Grid2d* u_next);
        float mitigation(float x, int borda);
        void applyAbsorptionBC();

        Domain d;
        Grid2d* u_current;
        Grid2d* u_next;
        Grid2d* f; 
        int frameRate;
};
