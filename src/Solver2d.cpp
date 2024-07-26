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
    // this->xs = int(this->X/2);
    // this->zs = int(this->Z/2);
    this-> xs = 500;
    this-> zs = 500;
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
    saveVTI(this->vel, this->outputFolder + "velocity_field.vti", "Velocity");
    if (rank == 0) {
        printf("Solving...\n");
    }
    auto start = chrono::high_resolution_clock::now();
    for (int k = 0; k <= Nt; k += 1){
        this->computeNext(k);
        // this->applyReynoldsBC(u_current, u_next);
        // this->applyAbsorptionBC();
        // this->computeNextInverted(k + 1);
        // this->applyReynoldsBC(u_next, u_current);
        // this->applyAbsorptionBC();
        if (k % frameRate == 0){
            string fileName = this->outputFolder + "output_data_" + to_string(rank) + "_" + to_string(k/frameRate) + ".vti";
            saveVTI(u_current, fileName, "Amplitude");
        }
        VecSwap(this->u_current, this->u_next);
    }
    auto final = chrono::high_resolution_clock::now();
    PetscFinalize();
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

void Solver2d::readInputVtkImageData(string input_file){
    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    printf("Rank: %d\n", rank);
    int* extent = new int[6];
    double* spacing = new double[3];
    vtkSmartPointer<vtkDataArray> data;
    if (rank == 0){
        vtkSmartPointer<vtkXMLImageDataReader> reader = vtkSmartPointer<vtkXMLImageDataReader>::New();
        reader->SetFileName(input_file.c_str());
        reader->Update();
        vtkSmartPointer<vtkImageData> imageData = reader->GetOutput();
        extent = imageData->GetExtent();
        spacing = imageData->GetSpacing();
        data = imageData->GetPointData()->GetScalars();
    }
    MPI_Bcast(extent, 6, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(spacing, 3, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
    this->Nx = extent[1] - extent[0] + 1;
    this->Nz = extent[3] - extent[2] + 1;
    this->dx = spacing[0];
    this->dz = spacing[1];
    DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DMDA_STENCIL_STAR, this->Nx, this->Nz, PETSC_DECIDE, PETSC_DECIDE, 1, 1, NULL, NULL, &this->da);
    DMDASetStencilWidth(this->da, STENCIL);
    DMSetFromOptions(this->da);
    DMSetUp(this->da);
    DMCreateGlobalVector(this->da, &this->vel);
    DMCreateLocalVector(this->da, &this->vel_local);
    int NxStart, NzStart, NxEnd, NzEnd;
    DMDAGetCorners(this->da, &NxStart, &NzStart, NULL, &NxEnd, &NzEnd, NULL);
    // VecGetArray(this->vel_local, &velArray);
    if (rank == 0){
        PetscScalar *velArray;
        velArray = new PetscScalar[Nx*Nz];
        PetscInt *indexArray;
        indexArray = new PetscInt[Nx*Nz];
        int k = 0;
        for (int i = NxStart; i < NxStart + NxEnd; i++){
            for (int j = NzStart; j < NzStart + NzEnd; j++){
                indexArray[k] = i*this->Nx + j;
                velArray[k] = (PetscScalar)data->GetTuple1(i*this->Nx + j);
                k++;
            }
        }
        VecSetValues(this->vel, this->Nz*this->Nx, indexArray, velArray, INSERT_VALUES);
    }
    VecAssemblyBegin(this->vel);
    VecAssemblyEnd(this->vel);
    DMGlobalToLocalBegin(this->da, this->vel, INSERT_VALUES, this->vel_local);
    DMGlobalToLocalEnd(this->da, this->vel, INSERT_VALUES, this->vel_local);
}

void Solver2d::saveVTI(Vec grid, string outputPath, string info){
    PetscInt NxStart, NzStart, NxEnd, NzEnd;
    DMDAGetCorners(da, &NxStart, &NzStart, NULL, &NxEnd, &NzEnd, NULL);
    vtkSmartPointer<vtkImageData> imageData = vtkSmartPointer<vtkImageData>::New();
    imageData->SetExtent(NxStart, NxStart + NxEnd - 1, NzStart, NzStart + NzEnd - 1, 0, 0);
    imageData->SetOrigin(NxStart*dx, NzStart*dz, 0);
    imageData->SetSpacing(dx, dz, 0);
    vtkSmartPointer<vtkFloatArray> data = vtkSmartPointer<vtkFloatArray>::New();
    data->SetNumberOfComponents(1);
    data->SetNumberOfTuples(Nx*Nz);
    data->SetName(info.c_str());
    PetscScalar *array;
    VecGetArray(grid, &array);
    for (int j = NzStart; j < NzStart + NzEnd; j++){
        for (int i = NxStart; i < NxStart + NxEnd; i++){
            data->SetTuple1(i + j*Nx, array[j*Nx + i]);
        }
    }
    VecRestoreArray(grid, &array);
    imageData->GetPointData()->SetScalars(data);
    vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();
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
    DMGlobalToLocalEnd(da,this->u_current,INSERT_VALUES,this->u_current_local);
    DMGlobalToLocalBegin(da,this->u_next,INSERT_VALUES,this->u_next_local);
    DMGlobalToLocalEnd(da,this->u_next,INSERT_VALUES,this->u_next_local);

    PetscScalar **u_currentArray, **u_nextArray, **velArray;
    DMDAVecGetArray(da, u_current_local, &u_currentArray);
    DMDAVecGetArray(da, u_next_local, &u_nextArray);
    DMDAVecGetArray(da, vel_local, &velArray);
    PetscInt NzStart, NxStart, NzEnd, NxEnd;
    DMDAGetCorners(da, &NxStart, &NzStart, NULL, &NxEnd, &NzEnd, NULL);
    NzEnd = NzStart + NzEnd;
    NxEnd = NxStart + NxEnd;
    if (NzStart < STENCIL){
        NzStart = STENCIL;
    }
    if (NxStart < STENCIL){
        NxStart = STENCIL;
    }
    if (NzEnd > this->Nz - STENCIL){
        NzEnd = NzEnd - STENCIL;
    }
    if (NxEnd > this->Nx - STENCIL){
        NxEnd = NxEnd - STENCIL;
    }
    for (int j = NzStart; j < NzEnd; j++){
        for (int i = NxStart; i < NxEnd; i++){
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

void Solver2d::computeNextInverted(int k){
    int sizez = Nz - STENCIL;
    int sizex = Nx - STENCIL;
    float val, courantNumber, const1, const2;

    DMGlobalToLocalBegin(da,this->u_current,INSERT_VALUES,this->u_current_local);
    DMGlobalToLocalEnd(da,this->u_current,INSERT_VALUES,this->u_current_local);
    DMGlobalToLocalBegin(da,this->u_next,INSERT_VALUES,this->u_next_local);
    DMGlobalToLocalEnd(da,this->u_next,INSERT_VALUES,this->u_next_local);

    PetscScalar **u_currentArray, **u_nextArray, **velArray;
    DMDAVecGetArray(da, u_current_local, &u_nextArray);
    DMDAVecGetArray(da, u_next_local, &u_currentArray);
    DMDAVecGetArray(da, vel_local, &velArray);
    PetscInt NzStart, NxStart, NzEnd, NxEnd;
    DMDAGetCorners(da, &NxStart, &NzStart, NULL, &NxEnd, &NzEnd, NULL);
    NzEnd = NzStart + NzEnd;
    NxEnd = NxStart + NxEnd;
    if (NzStart < STENCIL){
        NzStart = STENCIL;
    }
    if (NxStart < STENCIL){
        NxStart = STENCIL;
    }
    if (NzEnd > this->Nz - STENCIL){
        NzEnd = NzEnd - STENCIL;
    }
    if (NxEnd > this->Nx - STENCIL){
        NxEnd = NxEnd - STENCIL;
    }
    for (int j = NzStart; j < NzEnd; j++){
        for (int i = NxStart; i < NxEnd; i++){
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
    DMDAGetCorners(this->da, &NxStart, &NzStart, NULL, &Nx, &Nz, NULL);
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

