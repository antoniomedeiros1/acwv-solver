#include <iostream>
#include <fstream>
#include <cmath>
#include "math.h"
#include <string>
#include <chrono>

#include <vtkSmartPointer.h>
#include <vtkXMLImageDataReader.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkFloatArray.h>
#include <vtkDataArray.h>
#include <petscdm.h>
#include <petscdmda.h>

#define STENCIL 3

using namespace std;

class Solver2d{

    public:
        Solver2d(string input_file, string output_folder, int number_of_steps, float dt, int number_of_frames);
        ~Solver2d();
        void printParameters();
        void saveVTI(Vec u, string outputPath, string info);
        void saveVTIbin(Vec u, string outputPath, string info);
        void savePVTI(Vec u, string outputPath, string info);
        void solve();
    
    private: 
        void readInputTxt(string input_file);
        void readInputVtkImageData(string input_file);
        void computeNext(int k);
        void computeNextInverted(int k);
        float source(int x, int z, float k);
        void applyReynoldsBC(Vec u_current, Vec u_next);
        float mitigation(float x, int borda);
        void applyAbsorptionBC();

        PetscReal X, Y, Z, T, c, dx, dy, dz, dt;
        PetscInt Nx, Ny, Nz, Nt, xs, ys, zs;
        Vec vel;
        Vec vel_local;
        PetscReal cou, c1, c2;
        string outputFolder;
        DM da;
        Vec u_current;
        Vec u_current_local;
        Vec u_next;
        Vec u_next_local;
        int frameRate;
};
