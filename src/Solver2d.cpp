#include "../include/Solver2d.h"


Solver2d::Solver2d(string input_file, string output_folder, int number_of_steps, float dt, int number_of_frames){

    int nprocs;
    MPI_Comm_size(PETSC_COMM_WORLD, &nprocs);
    this->outputFolder = output_folder;
    printf("Reading input file: %s\n", input_file.c_str());
    this->readInputVtkImageData(input_file);
    this->dt = dt;
    this->Nt = number_of_steps;
    this->X = this->Nx * this->dx;
    this->Z = this->Nz * this->dz;
    this->T = this->Nt * this->dt;
    this->frameRate = this->Nt/number_of_frames;
    this->xs = int(this->X/2);
    this->zs = int(this->Z/2);
    DMCreateGlobalVector(this->da, &this->u_current);
    DMCreateLocalVector(this->da, &this->u_current_local);
    DMCreateGlobalVector(this->da, &this->u_next);
    DMCreateLocalVector(this->da, &this->u_next_local);
    VecSet(this->u_current, 0.0);
    VecAssemblyBegin(this->u_current);
    VecAssemblyEnd(this->u_current);
    VecSet(this->u_next, 0.0);
    VecAssemblyBegin(this->u_next);
    VecAssemblyEnd(this->u_next);
}

Solver2d::~Solver2d(){}

void Solver2d::solve(){
    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    if (rank == 0) {
        printParameters();
        printf("Saving velocity field...\n");
    }
    // saveVTI(this->vel, this->outputFolder + "velocity_field.vti", "Velocity");
    if (rank == 0) {
        printf("Solving...\n");
    }
    auto start = chrono::high_resolution_clock::now();
    for (int k = 0; k <= Nt; k += 2){
        this->computeNext(k);
        // this->applyReynoldsBC(u_current, u_next);
        // this->applyAbsorptionBC();
        this->computeNextInverted(k + 1);
        // this->applyReynoldsBC(u_next, u_current);
        // this->applyAbsorptionBC();
        if (k % frameRate == 0){
            string fileName = this->outputFolder + "output_data_" + to_string(rank) + "_" + to_string(k/frameRate) + ".vti";
            // saveVTI(u_current, fileName, "Amplitude");
        }
    }
    auto final = chrono::high_resolution_clock::now();
    printf("File saved\n");
    chrono::duration<double> interval = final - start;
    printf("Elapsed time: %f seconds\n", interval.count());
}

void Solver2d::printParameters(){
    printf("\nSimulation paremeters:\n");
    cout << "X = " << this->X << "m\n";
    cout << "Z = " << this->Z << "m\n"; 
    cout << "T = " << this->T << "s\n";
    cout << "dx = " << this->dx << "m\n";
    cout << "dt = " << this->dt << "s\n";
}

// void Solver2d::readInputTxt(string input_file){
//     ifstream myfile;
//     PetscScalar **velArray;
//     myfile.open(input_file);
//     if(myfile.is_open()){
//         myfile >> this->d.Nx; 
//         myfile >> this->d.Nz;
//         myfile >> this->d.Nt;
//         myfile >> this->d.dx;
//         this->d.dz = this->d.dx;
//         myfile >> this->d.dt;
//         DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DMDA_STENCIL_STAR, this->d.Nz, this->d.Nx, PETSC_DECIDE, PETSC_DECIDE, 1, 1, NULL, NULL, &this->da);
//         DMDASetStencilWidth(this->da, STENCIL);
//         DMCreateGlobalVector(this->da, &this->d.vel);
//         DMDAVecGetArray(this->da, this->d.vel, &velArray);
//         for (PetscInt j = 0; j < this->d.Nz; j++){
//             for (PetscInt i = 0; i < this->d.Nx; i++){
//                 PetscReal val;
//                 myfile >> val;
//                 velArray[j][i] = val;
//             }
//         }
//         DMDAVecRestoreArray(this->da, this->d.vel, &velArray);
//         myfile.close();
//     } else {
//         printf("Failed to open file\n");
//         exit(1);
//     }
// }

void Solver2d::readInputVtkImageData(string input_file){
    int rank;
    int size[2];
    double dxz[2];
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    PetscScalar *velArray;
    PetscInt *indexArray;
    if (rank == 0){
        vtkSmartPointer<vtkDataArray> data;
        vtkSmartPointer<vtkXMLImageDataReader> reader;
        vtkSmartPointer<vtkImageData> imageData;
        reader = vtkSmartPointer<vtkXMLImageDataReader>::New();
        reader->SetFileName(input_file.c_str());
        reader->Update();
        imageData = reader->GetOutput();
        int* extent = imageData->GetExtent();
        double* spacing = imageData->GetSpacing();
        size[0] = extent[1] - extent[0] + 1;
        size[1] = extent[3] - extent[2] + 1;
        dxz[0] = spacing[0];
        dxz[1] = spacing[1];
        data = imageData->GetPointData()->GetScalars();
        velArray = new PetscScalar[size[0]*size[1]];
        indexArray = new PetscInt[size[0]*size[1]];
        int k = 0;
        for (int i = 0; i < this->Nz; i++){
            for (int j = 0; j < this->Nx; j++){
                indexArray[k] = i*this->Nx + j;
                velArray[k] = (PetscScalar)data->GetTuple1(i*this->Nx + j);
                k++;
            }
        }
    }
    MPI_Bcast(size, 2, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(dxz, 2, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
    this->Nx = size[0];
    this->Nz = size[1];
    this->dx = dxz[0];
    this->dz = dxz[1];
    // PetscScalar **velArray;
    DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DMDA_STENCIL_STAR, this->Nz, this->Nx, PETSC_DECIDE, PETSC_DECIDE, 1, 1, NULL, NULL, &this->da);
    printf("Rank: %d\n", rank);
    DMDASetStencilWidth(this->da, 2);
    DMSetUp(this->da);
    DMCreateGlobalVector(this->da, &this->vel);
    DMCreateLocalVector(this->da, &this->vel_local);
    if (rank == 0){
        // PetscScalar velArray[this->Nz * this->Nx];
        // PetscInt indexArray[this->Nz * this->Nx];
        VecSetValues(this->vel, this->Nz*this->Nx, indexArray, velArray, INSERT_VALUES);
        delete[] velArray;
        delete[] indexArray;
    }
    VecAssemblyBegin(this->vel);
    VecAssemblyEnd(this->vel);
    DMGlobalToLocalBegin(this->da, this->vel, INSERT_VALUES, this->vel_local);
    DMGlobalToLocalEnd(this->da, this->vel, INSERT_VALUES, this->vel_local);
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

void Solver2d::saveVTI(Vec grid, string outputPath, string info){
    vtkSmartPointer<vtkImageData> imageData = vtkSmartPointer<vtkImageData>::New();
    PetscInt NxStart, NzStart, Nx, Nz;
    DMDAGetCorners(da, &NzStart, &NxStart, NULL, &Nz, &Nx, NULL);
    imageData->SetExtent(0, Nx - 1, 0, Nz - 1, 0, 0);
    imageData->SetOrigin(0*dx, (Nz - 1)*dz, 0);
    imageData->SetSpacing(dx, dz, 0);
    vtkSmartPointer<vtkFloatArray> data = vtkSmartPointer<vtkFloatArray>::New();
    data->SetNumberOfComponents(1);
    data->SetNumberOfTuples(Nx*Nz);
    data->SetName(info.c_str());
    PetscScalar *array;
    VecGetArray(grid, &array);
    for (int j = NzStart; j < NzStart + Nz; j++){
        for (int i = NxStart; i < NxStart + Nx; i++){
            data->SetTuple1(i + j*Nx, array[j*Nx + i]);
        }
    }
    VecRestoreArray(grid, &array);
    imageData->GetPointData()->SetScalars(data);
    vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();
    // ascii
    writer->SetDataModeToAscii();
    writer->SetFileName(outputPath.c_str());
    writer->SetInputData(imageData);
    writer->Write();
}

void Solver2d::saveVTIbin(Vec grid, string outputPath, string info){
    vtkSmartPointer<vtkImageData> imageData = vtkSmartPointer<vtkImageData>::New();
    imageData->SetExtent(STENCIL, Nx - 1 - STENCIL, STENCIL, Nz - 1 - STENCIL, 0, 0);
    imageData->SetOrigin(STENCIL*dx, (Nz - 1)*dz, 0);
    imageData->SetSpacing(dx, dz, 0);
    vtkSmartPointer<vtkFloatArray> data = vtkSmartPointer<vtkFloatArray>::New();
    data->SetNumberOfComponents(1);
    data->SetNumberOfTuples(Nx*Nz);
    data->SetName(info.c_str());
    PetscScalar *array;
    VecGetArray(grid, &array);
    for (int j = STENCIL; j < Nz - STENCIL; j++){
        for (int i = STENCIL; i < Nx - STENCIL; i++){
            data->SetTuple1(i + j*Nx, array[j*Nx + i]);
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
    if (x != (int)(this->xs/this->dx) || z!= (int)(this->zs/this->dz) || k*this->dt > 0.5){
        return 0;
    } 
    float td = k*this->dt - ((2.0f*sqrtf(M_PI))/fcorte);  
    float fc = (fcorte/(3.0f*sqrtf(M_PI)));
    return (1.0f - 2.0f * M_PI * powf(M_PI * fc * td, 2.0f))/powf(M_E, M_PI*powf((M_PI*fc*td), 2.0f));
}

void Solver2d::computeNext(int k){
    int sizez = Nz - STENCIL;
    int sizex = Nx - STENCIL;
    float val, courantNumber, const1, const2;

    DMGlobalToLocalBegin(da,this->u_current,INSERT_VALUES,this->u_current_local);
    DMGlobalToLocalEnd(da,this->u_next,INSERT_VALUES,this->u_next_local);

    PetscScalar **u_currentArray, **u_nextArray, **velArray;
    DMDAVecGetArray(da, u_current_local, &u_currentArray);
    DMDAVecGetArray(da, u_next_local, &u_nextArray);
    DMDAVecGetArray(da, vel_local, &velArray);
    PetscInt NzStart, NxStart, NzSize, NxSize;
    DMDAGetCorners(da, &NzStart, &NxStart, NULL, &NzSize, &NxSize, NULL);
    for (int j = NzStart; j < NzStart + NzSize; j++){
        for (int i = NxStart; i < NxStart + NxSize; i++){
            if (j < 2 || j >= Nz - 2 || i < 2 || i >= Nx - 2){
                u_currentArray[j][i] = 0;
                u_nextArray[j][i] = 0;
                continue;
            }
            courantNumber = dt*velArray[j][i]/dx;
            const1 = (powf(courantNumber, 2.0f)/12.0f);
            const2 = powf(velArray[j][i]*dt, 2.0f);
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
    DMDAVecRestoreArray(da, u_current_local, &u_nextArray);
    DMDAVecRestoreArray(da, u_next_local, &u_currentArray);

    DMLocalToGlobalBegin(da,u_current_local,INSERT_VALUES,u_current);
    DMLocalToGlobalEnd(da,u_current_local,INSERT_VALUES,u_current);
    DMLocalToGlobalBegin(da,u_next_local,INSERT_VALUES,u_next);
    DMLocalToGlobalEnd(da,u_next_local,INSERT_VALUES,u_next);
}

void Solver2d::computeNextInverted(int k){
    int sizez = Nz - STENCIL;
    int sizex = Nx - STENCIL;
    float val, courantNumber, const1, const2;

    DMGlobalToLocalBegin(da,this->u_current,INSERT_VALUES,this->u_current_local);
    DMGlobalToLocalEnd(da,this->u_next,INSERT_VALUES,this->u_next_local);

    PetscScalar **u_currentArray, **u_nextArray, **velArray;
    DMDAVecGetArray(da, u_current_local, &u_nextArray);
    DMDAVecGetArray(da, u_next_local, &u_currentArray);
    DMDAVecGetArray(da, vel_local, &velArray);
    PetscInt NzStart, NxStart, NzSize, NxSize;
    DMDAGetCorners(da, &NzStart, &NxStart, NULL, &NzSize, &NxSize, NULL);
    for (int j = NzStart; j < NzStart + NzSize; j++){
        for (int i = NxStart; i < NxStart + NxSize; i++){
            if (j < 2 || j >= Nz - 2 || i < 2 || i >= Nx - 2){
                u_currentArray[j][i] = 0;
                u_nextArray[j][i] = 0;
                continue;
            }
            courantNumber = dt*velArray[j][i]/dx;
            const1 = (powf(courantNumber, 2.0f)/12.0f);
            const2 = powf(velArray[j][i]*dt, 2.0f);
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
    DMDAVecRestoreArray(da, u_current_local, &u_currentArray);
    DMDAVecRestoreArray(da, u_next_local, &u_nextArray);

    DMLocalToGlobalBegin(da,u_current_local,INSERT_VALUES,u_current);
    DMLocalToGlobalEnd(da,u_current_local,INSERT_VALUES,u_current);
    DMLocalToGlobalBegin(da,u_next_local,INSERT_VALUES,u_next);
    DMLocalToGlobalEnd(da,u_next_local,INSERT_VALUES,u_next);
}

void Solver2d::applyReynoldsBC(Vec u_current, Vec u_next){
    PetscScalar **u_currentArray, **u_nextArray, **velArray;
    DMDAVecGetArray(this->da, u_current, &u_currentArray);
    DMDAVecGetArray(this->da, u_next, &u_nextArray);
    DMDAVecGetArray(this->da, this->vel, &velArray);
    PetscInt NzStart, NxStart, Nz, Nx;
    DMDAGetCorners(this->da, &NzStart, &NxStart, NULL, &Nz, &Nx, NULL);
    for(int j = STENCIL; j < Nz - STENCIL; j++){
        for(int i = STENCIL; i <= STENCIL + 1; i++){
            // float courantNumber = dt * vel->get(j, i)/dx;
            // (*u_next)(j, i) = (*u_current)(j, i) + courantNumber*((*u_current)(j, i + 1) - (*u_current)(j, i));
            float courantNumber = dt * velArray[j][i]/dx;
            u_nextArray[j][i] = u_currentArray[j][i] + courantNumber*(u_currentArray[j][i + 1] - u_currentArray[j][i]);
        }
    }
    for(int j = STENCIL; j < Nz - STENCIL; j++){
        for(int i = Nx - STENCIL - 1; i <= Nx - STENCIL; i++){
            // float courantNumber = dt * vel->get(j, i)/dx;
            // (*u_next)(j, i) = (*u_current)(j, i) - courantNumber*((*u_current)(j, i) - (*u_current)(j, i - 1));
            float courantNumber = dt * velArray[j][i]/dx;
            u_nextArray[j][i] = u_currentArray[j][i] - courantNumber*(u_currentArray[j][i] - u_currentArray[j][i - 1]);
        }
    }
    for(int i = STENCIL; i < Nx - STENCIL; i++){
        for(int j = STENCIL; j <= STENCIL + 1; j++){
            // float courantNumber = dt * vel->get(j, i)/dx;
            // (*u_next)(j, i) = (*u_current)(j, i) + courantNumber*((*u_current)(j + 1, i) - (*u_current)(j, i));
            float courantNumber = dt * velArray[j][i]/dx;
            u_nextArray[j][i] = u_currentArray[j][i] + courantNumber*(u_currentArray[j + 1][i] - u_currentArray[j][i]);
        }
    }
    for(int i = STENCIL; i < Nx - STENCIL; i++){
        for(int j = Nz - STENCIL - 1; j <= Nz - STENCIL; j++){
            // float courantNumber = dt * vel->get(j, i)/dx;
            // (*u_next)(j, i) = (*u_current)(j, i) - courantNumber*((*u_current)(j, i) - (*u_current)(j - 1, i));
            float courantNumber = dt * velArray[j][i]/dx;
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
        for(int i = STENCIL; i <= Nx - STENCIL; i++){
            // u_current->set(j, i, u_current->get(j, i)*mitigation(j, border));
            // u_next->set(j, i, u_next->get(j, i)*mitigation(j, border));
            u_currentArray[j][i] = u_currentArray[j][i]*mitigation(j, border);
            u_nextArray[j][i] = u_nextArray[j][i]*mitigation(j, border);
        }
    }
    for(int j = Nz - border; j <= Nz - STENCIL; j++){
        for(int i = STENCIL; i <= Nx - STENCIL; i++){
            // u_current->set(j, i, u_current->get(j, i)*mitigation(Nz - j, border));
            // u_next->set(j, i, u_next->get(j, i)*mitigation(Nz - j, border));
            u_currentArray[j][i] = u_currentArray[j][i]*mitigation(Nz - j, border);
            u_nextArray[j][i] = u_nextArray[j][i]*mitigation(Nz - j, border);
        }
    }
    for(int j = border; j < Nz - border; j++){
        for(int i = STENCIL; i <= border; i++){
            // u_current->set(j,i, u_current->get(j, i)*mitigation(i, border));
            // u_next->set(j,i, u_next->get(j, i)*mitigation(i, border));
            u_currentArray[j][i] = u_currentArray[j][i]*mitigation(i, border);
            u_nextArray[j][i] = u_nextArray[j][i]*mitigation(i, border);
        }
    }
    for(int j = border; j < Nz - border; j++){
        for(int i = Nx - border; i <= Nx - STENCIL; i++){
            // u_current->set(j,i, u_current->get(j, i)*mitigation((Nx - i), border));
            // u_next->set(j,i, u_next->get(j, i)*mitigation((Nx - i), border));
            u_currentArray[j][i] = u_currentArray[j][i]*mitigation((Nx - i), border);
            u_nextArray[j][i] = u_nextArray[j][i]*mitigation((Nx - i), border);
        }
    }
}

