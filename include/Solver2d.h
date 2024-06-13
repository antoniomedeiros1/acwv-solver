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

#include "Domain.h"

#define STENCIL 3

using namespace std;

class Solver2d{

    public:
        Solver2d(string input_file, string output_folder, int number_of_steps, float dt, int number_of_frames);
        ~Solver2d();
        void printParameters();
        void saveVTI(Domain d, Vec u, string outputPath, string info);
        void saveVTIbin(Domain d, Vec u, string outputPath, string info);
        void savePVTI(Domain d, Vec u, string outputPath, string info);
        void solve();
    
    private: 
        void readInputTxt(string input_file);
        void readInputVtkImageData(string input_file);
        void computeNext(Domain d, int k);
        void computeNextInverted(Domain d, int k);
        float source(int x, int z, float k);
        void applyReynoldsBC(Vec u_current, Vec u_next);
        float mitigation(float x, int borda);
        void applyAbsorptionBC();

        string outputFolder;
        Domain d;
        DM da;
        Vec u_current;
        Vec u_current_local;
        Vec u_next;
        Vec u_next_local;
        int frameRate;
};
