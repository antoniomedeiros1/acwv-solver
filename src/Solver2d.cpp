#include "../include/Solver2d.h"


Solver2d::Solver2d(string input_file, string output_folder, int number_of_steps, float dt, int number_of_frames){
    PetscInitialize(NULL, NULL, NULL, NULL);
    printf("Reading input file: %s\n", input_file.c_str());
    this->outputFolder = output_folder;
    this->d = Domain();
    this->d.dt = dt;
    this->d.Nt = number_of_steps;
    this->readInputVtkImageData(input_file);
    this->d.X = this->d.Nx * this->d.dx;
    this->d.Z = this->d.Nz * this->d.dz;
    this->d.T = this->d.Nt * this->d.dt;
    this->frameRate = this->d.Nt/number_of_frames;
    this->d.xs = int(this->d.X/2);
    this->d.zs = int(this->d.Z/2);
    DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DMDA_STENCIL_STAR, this->d.Nz, this->d.Nx, PETSC_DECIDE, PETSC_DECIDE, 1, 1, NULL, NULL, &this->da);
    DMCreateGlobalVector(this->da, &this->u_current);
    DMCreateGlobalVector(this->da, &this->u_next);
    VecSet(this->u_current, 0.0);
    VecAssemblyBegin(this->u_current);
    VecAssemblyEnd(this->u_current);
    VecSet(this->u_next, 0.0);
    VecAssemblyBegin(this->u_next);
    VecAssemblyEnd(this->u_next);
}

Solver2d::~Solver2d(){}

void Solver2d::solve(){
    printParameters();
    printf("Saving velocity field...\n");
    saveVTI(this->d, this->d.vel, this->outputFolder + "velocity_field", "Velocity");
    printf("Solving...\n");
    auto start = chrono::high_resolution_clock::now();
    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    for (int k = 0; k <= d.Nt; k += 2){
        this->computeNext(d, u_current, u_next, k);
        this->applyReynoldsBC(u_current, u_next);
        this->applyAbsorptionBC();
        this->computeNext(d, u_next, u_current, k + 1);
        this->applyReynoldsBC(u_next, u_current);
        this->applyAbsorptionBC();
        if (k % frameRate == 0 && rank == 0){
            string fileName = this->outputFolder + "output_data" + to_string(k/frameRate) + ".vti";
            saveVTI(d, u_current, fileName, "Amplitude");
        }
    }
    auto final = chrono::high_resolution_clock::now();
    PetscFinalize();
    printf("File saved\n");
    chrono::duration<double> interval = final - start;
    printf("Elapsed time: %f seconds\n", interval.count());
}

void Solver2d::printParameters(){
    printf("\nSimulation paremeters:\n");
    cout << "X = " << this->d.X << "m\n";
    cout << "Z = " << this->d.Z << "m\n"; 
    cout << "T = " << this->d.T << "s\n";
    cout << "dx = " << this->d.dx << "m\n";
    cout << "dt = " << this->d.dt << "s\n";
}

void Solver2d::readInputTxt(string input_file){
    ifstream myfile;
    PetscScalar **velArray;
    myfile.open(input_file);
    if(myfile.is_open()){
        myfile >> this->d.Nx; 
        myfile >> this->d.Nz;
        myfile >> this->d.Nt;
        myfile >> this->d.dx;
        this->d.dz = this->d.dx;
        myfile >> this->d.dt;
        DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DMDA_STENCIL_STAR, this->d.Nz, this->d.Nx, PETSC_DECIDE, PETSC_DECIDE, 1, 1, NULL, NULL, &this->d.da);
        DMCreateGlobalVector(this->d.da, &this->d.vel);
        DMDAVecGetArray(this->d.da, this->d.vel, &velArray);
        for (PetscInt j = 0; j < this->d.Nz; j++){
            for (PetscInt i = 0; i < this->d.Nx; i++){
                PetscReal val;
                myfile >> val;
                velArray[j][i] = val;
            }
        }
        DMDAVecRestoreArray(this->d.da, this->d.vel, &velArray);
        myfile.close();
    } else {
        printf("Failed to open file\n");
        exit(1);
    }
}

void Solver2d::readInputVtkImageData(string input_file){
    vtkSmartPointer<vtkXMLImageDataReader> reader = vtkSmartPointer<vtkXMLImageDataReader>::New();
    reader->SetFileName(input_file.c_str());
    reader->Update();
    vtkSmartPointer<vtkImageData> imageData = reader->GetOutput();
    int* extent = imageData->GetExtent();
    double* spacing = imageData->GetSpacing();
    this->d.Nx = extent[1] - extent[0] + 1;
    this->d.Nz = extent[3] - extent[2] + 1;
    this->d.dx = spacing[0];
    this->d.dz = spacing[1];
    vtkSmartPointer<vtkDataArray> data = imageData->GetPointData()->GetScalars();
    PetscScalar **velArray;
    DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DMDA_STENCIL_STAR, this->d.Nz, this->d.Nx, PETSC_DECIDE, PETSC_DECIDE, 1, 1, NULL, NULL, &this->d.da);
    DMCreateGlobalVector(this->d.da, &this->d.vel);
    DMDAVecGetArray(this->d.da, this->d.vel, &velArray);
    for (int j = 0; j < this->d.Nz; j++){
        for (int i = 0; i < this->d.Nx; i++){
            PetscReal val = data->GetTuple1(i + j*this->d.Nx);
            velArray[j][i] = val;
        }
    }
    DMDAVecRestoreArray(this->d.da, this->d.vel, &velArray);
}

// void Solver2d::saveVTI(Domain d, Vec grid, string outputPath, string info){
//     ofstream myfile;
//     myfile.open(outputPath + ".vti");
//     if(myfile.is_open()){
//         myfile << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
//         myfile << "  <ImageData WholeExtent= \"" <<  STENCIL << " " << d.Nx - 1 - STENCIL << " " << STENCIL << " " << d.Nz - 1 - STENCIL << " " << 0 << " " << 0 << "\" ";
//         myfile << "Origin = \"" << STENCIL << " " << d.Nz - 1 << " " << 0 << "\" ";
//         myfile << "Spacing = \"" << d.dx << " " << d.dz << " " << 0 << "\">\n";
//         myfile << "    <Piece Extent = \"" << STENCIL << " " << d.Nx - 1 - STENCIL << " " << STENCIL << " " << d.Nz - 1 - STENCIL << " " << 0 << " " << 0 << "\">\n";
//         myfile << "      <PointData Scalars=\"" + info + "\">\n";
//         myfile << "        <DataArray type=\"Float32\" Name=\"" + info + "\" format=\"ascii\">\n";
//         for (int j = STENCIL; j < d.Nz - STENCIL; j++){
//             for (int i = STENCIL; i < d.Nx - STENCIL; i++){
//                 PetscReal val;
//                 VecGetValues(grid, 1, &j*d.Nx + i, &val);
//                 myfile << val << " ";
//             }
//         }
//         myfile << "\n        </DataArray>";
//         myfile << "\n      </PointData>";
//         myfile << "\n    </Piece>";
//         myfile << "\n  </ImageData>";
//         myfile << "\n</VTKFile>";
//     } else {
//         cout << "Failed to create " << outputPath << ".vti" << endl;
//     }
// }

void Solver2d::saveVTI(Domain d, Vec grid, string outputPath, string info){
    vtkSmartPointer<vtkImageData> imageData = vtkSmartPointer<vtkImageData>::New();
    imageData->SetExtent(STENCIL, d.Nx - 1 - STENCIL, STENCIL, d.Nz - 1 - STENCIL, 0, 0);
    imageData->SetOrigin(STENCIL*d.dx, (d.Nz - 1)*d.dz, 0);
    imageData->SetSpacing(d.dx, d.dz, 0);
    vtkSmartPointer<vtkFloatArray> data = vtkSmartPointer<vtkFloatArray>::New();
    data->SetNumberOfComponents(1);
    data->SetNumberOfTuples(d.Nx*d.Nz);
    data->SetName(info.c_str());
    PetscScalar *array;
    VecGetArray(grid, &array);
    for (int j = STENCIL; j < d.Nz - STENCIL; j++){
        for (int i = STENCIL; i < d.Nx - STENCIL; i++){
            data->SetTuple1(i + j*d.Nx, array[j*d.Nx + i]);
        }
    }
    VecRestoreArray(grid, &array);
    imageData->GetPointData()->SetScalars(data);
    vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();
    writer->SetFileName(outputPath.c_str());
    writer->SetInputData(imageData);
    writer->Write();
}

void Solver2d::saveVTIbin(Domain d, Vec grid, string outputPath, string info){
    vtkSmartPointer<vtkImageData> imageData = vtkSmartPointer<vtkImageData>::New();
    imageData->SetExtent(STENCIL, d.Nx - 1 - STENCIL, STENCIL, d.Nz - 1 - STENCIL, 0, 0);
    imageData->SetOrigin(STENCIL*d.dx, (d.Nz - 1)*d.dz, 0);
    imageData->SetSpacing(d.dx, d.dz, 0);
    vtkSmartPointer<vtkFloatArray> data = vtkSmartPointer<vtkFloatArray>::New();
    data->SetNumberOfComponents(1);
    data->SetNumberOfTuples(d.Nx*d.Nz);
    data->SetName(info.c_str());
    PetscScalar *array;
    VecGetArray(grid, &array);
    for (int j = STENCIL; j < d.Nz - STENCIL; j++){
        for (int i = STENCIL; i < d.Nx - STENCIL; i++){
            data->SetTuple1(i + j*d.Nx, array[j*d.Nx + i]);
        }
    }
    VecRestoreArray(grid, &array);
    imageData->GetPointData()->SetScalars(data);
    vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();
    writer->SetFileName(outputPath.c_str());
    writer->SetInputData(imageData);
    writer->SetDataModeToBinary();
    writer->Write();
}

float Solver2d::source(int x, int z, float k){
    float fcorte = 40;
    if (x != (int)(this->d.xs/this->d.dx) || z!= (int)(this->d.zs/this->d.dz) || k*this->d.dt > 0.5){
        return 0;
    } 
    float td = k*this->d.dt - ((2.0f*sqrtf(M_PI))/fcorte);  
    float fc = (fcorte/(3.0f*sqrtf(M_PI)));
    return (1.0f - 2.0f * M_PI * powf(M_PI * fc * td, 2.0f))/powf(M_E, M_PI*powf((M_PI*fc*td), 2.0f));
}

void Solver2d::computeNext(Domain d, Vec u_current, Vec u_next, int k){
    int sizez = d.Nz - STENCIL;
    int sizex = d.Nx - STENCIL;
    float val, courantNumber, const1, const2;
    PetscScalar **u_currentArray, **u_nextArray, **velArray;
    DMDAVecGetArray(d.da, u_current, &u_currentArray);
    DMDAVecGetArray(d.da, u_next, &u_nextArray);
    DMDAVecGetArray(d.da, d.vel, &velArray);
    for (int j = STENCIL; j <= sizez; j++){
        for (int i = STENCIL; i <= sizex; i++){
            courantNumber = d.dt*velArray[j][i]/d.dx;
            const1 = (powf(courantNumber, 2.0f)/12.0f);
            const2 = powf(velArray[j][i]*d.dt, 2.0f);
            val = 
            const1 *
            (
                -1*(u_currentArray[j][i - 2] + u_currentArray[j - 2][i]) + 
                16*(u_currentArray[j][i - 1] + u_currentArray[j - 1][i]) - 
                60* u_currentArray[j][i] +
                16*(u_currentArray[j][i + 1] + u_currentArray[j + 1][i]) -
                   (u_currentArray[j][i + 2] + u_currentArray[j + 2][i]) 
            ) 
            + 2*u_currentArray[j][i] - u_nextArray[j][i] - const2 * this->source(i, j, k);
            u_nextArray[j][i] = val;
        }
    }
    DMDAVecRestoreArray(d.da, u_current, &u_currentArray);
    DMDAVecRestoreArray(d.da, u_next, &u_nextArray);
}

void Solver2d::applyReynoldsBC(Vec u_current, Vec u_next){
    PetscScalar **u_currentArray, **u_nextArray, **velArray;
    DMDAVecGetArray(this->da, u_current, &u_currentArray);
    DMDAVecGetArray(this->da, u_next, &u_nextArray);
    DMDAVecGetArray(this->da, this->d.vel, &velArray);
    for(int j = STENCIL; j < d.Nz - STENCIL; j++){
        for(int i = STENCIL; i <= STENCIL + 1; i++){
            // float courantNumber = d.dt * d.vel->get(j, i)/d.dx;
            // (*u_next)(j, i) = (*u_current)(j, i) + courantNumber*((*u_current)(j, i + 1) - (*u_current)(j, i));
            float courantNumber = d.dt * velArray[j][i]/d.dx;
            u_nextArray[j][i] = u_currentArray[j][i] + courantNumber*(u_currentArray[j][i + 1] - u_currentArray[j][i]);
        }
    }
    for(int j = STENCIL; j < d.Nz - STENCIL; j++){
        for(int i = d.Nx - STENCIL - 1; i <= d.Nx - STENCIL; i++){
            // float courantNumber = d.dt * d.vel->get(j, i)/d.dx;
            // (*u_next)(j, i) = (*u_current)(j, i) - courantNumber*((*u_current)(j, i) - (*u_current)(j, i - 1));
            float courantNumber = d.dt * velArray[j][i]/d.dx;
            u_nextArray[j][i] = u_currentArray[j][i] - courantNumber*(u_currentArray[j][i] - u_currentArray[j][i - 1]);
        }
    }
    for(int i = STENCIL; i < d.Nx - STENCIL; i++){
        for(int j = STENCIL; j <= STENCIL + 1; j++){
            // float courantNumber = d.dt * d.vel->get(j, i)/d.dx;
            // (*u_next)(j, i) = (*u_current)(j, i) + courantNumber*((*u_current)(j + 1, i) - (*u_current)(j, i));
            float courantNumber = d.dt * velArray[j][i]/d.dx;
            u_nextArray[j][i] = u_currentArray[j][i] + courantNumber*(u_currentArray[j + 1][i] - u_currentArray[j][i]);
        }
    }
    for(int i = STENCIL; i < d.Nx - STENCIL; i++){
        for(int j = d.Nz - STENCIL - 1; j <= d.Nz - STENCIL; j++){
            // float courantNumber = d.dt * d.vel->get(j, i)/d.dx;
            // (*u_next)(j, i) = (*u_current)(j, i) - courantNumber*((*u_current)(j, i) - (*u_current)(j - 1, i));
            float courantNumber = d.dt * velArray[j][i]/d.dx;
            u_nextArray[j][i] = u_currentArray[j][i] - courantNumber*(u_currentArray[j][i] - u_currentArray[j - 1][i]);
        }
    }
    DMDAVecRestoreArray(this->da, u_current, &u_currentArray);
    DMDAVecRestoreArray(this->da, u_next, &u_nextArray);
}

float Solver2d::mitigation(float x, int border){
    float fat = 0.0055f;
    return expf(-(powf(fat*(border - x), 2.0f)));
}

void Solver2d::applyAbsorptionBC(){
    int border = 25;
    PetscScalar **u_currentArray, **u_nextArray;
    DMDAVecGetArray(this->da, this->u_current, &u_currentArray);
    DMDAVecGetArray(this->da, this->u_next, &u_nextArray);
    for(int j = STENCIL; j < border; j++){
        for(int i = STENCIL; i <= d.Nx - STENCIL; i++){
            // u_current->set(j, i, u_current->get(j, i)*mitigation(j, border));
            // u_next->set(j, i, u_next->get(j, i)*mitigation(j, border));
            u_currentArray[j][i] = u_currentArray[j][i]*mitigation(j, border);
            u_nextArray[j][i] = u_nextArray[j][i]*mitigation(j, border);
        }
    }
    for(int j = d.Nz - border; j <= d.Nz - STENCIL; j++){
        for(int i = STENCIL; i <= d.Nx - STENCIL; i++){
            // u_current->set(j, i, u_current->get(j, i)*mitigation(d.Nz - j, border));
            // u_next->set(j, i, u_next->get(j, i)*mitigation(d.Nz - j, border));
            u_currentArray[j][i] = u_currentArray[j][i]*mitigation(d.Nz - j, border);
            u_nextArray[j][i] = u_nextArray[j][i]*mitigation(d.Nz - j, border);
        }
    }
    for(int j = border; j < d.Nz - border; j++){
        for(int i = STENCIL; i <= border; i++){
            // u_current->set(j,i, u_current->get(j, i)*mitigation(i, border));
            // u_next->set(j,i, u_next->get(j, i)*mitigation(i, border));
            u_currentArray[j][i] = u_currentArray[j][i]*mitigation(i, border);
            u_nextArray[j][i] = u_nextArray[j][i]*mitigation(i, border);
        }
    }
    for(int j = border; j < d.Nz - border; j++){
        for(int i = d.Nx - border; i <= d.Nx - STENCIL; i++){
            // u_current->set(j,i, u_current->get(j, i)*mitigation((d.Nx - i), border));
            // u_next->set(j,i, u_next->get(j, i)*mitigation((d.Nx - i), border));
            u_currentArray[j][i] = u_currentArray[j][i]*mitigation((d.Nx - i), border);
            u_nextArray[j][i] = u_nextArray[j][i]*mitigation((d.Nx - i), border);
        }
    }
}

