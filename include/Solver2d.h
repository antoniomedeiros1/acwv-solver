#include <iostream>
#include <fstream>
#include <cmath>
#include "math.h"
#include <string>
#include <chrono>

#include <vtkSmartPointer.h>
#include <vtkXMLImageDataReader.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkFloatArray.h>
#include <vtkDataArray.h>

#include "Domain.h"

#define STENCIL 3

using namespace std;

class Solver2d{

    public:
        Solver2d(string input_file, string output_folder, int number_of_steps, float dt, int number_of_frames);
        ~Solver2d();
        void printParameters();
        void saveVTI(Domain d, Grid2d* u, string outputPath, string info);
        void saveVTIbin(Domain d, Grid2d* u, string outputPath, string info);
        void solve();
    
    private: 
        void readInputTxt(string input_file);
        void readInputVtkImageData(string input_file);
        void computeNext(Domain d, Grid2d* u_current, Grid2d* u_next, int k);
        void computeNext(Domain d, Grid2d* u_current, Grid2d* u_next);
        float source(int x, int z, float k);
        void applyReynoldsBC(Grid2d* u_current, Grid2d* u_next);
        float mitigation(float x, int borda);
        void applyAbsorptionBC();

        string outputFolder;
        Domain d;
        Grid2d* u_current;
        Grid2d* u_next;
        Grid2d* f; 
        int frameRate;
};
