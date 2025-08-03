#include "../include/Solver2d.h"


Solver2d::Solver2d(
        string input_file,
        string output_folder,
        int xs,
        int zs,
        int number_of_steps,
        float dt,
        int number_of_frames
){
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
    this->xs = xs;
    this->zs = zs;
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
    writeVTI(this->vel, this->outputFolder + "velocity_field_" + to_string(rank) + ".vti", "Velocity");
    if (rank == 0) {
        printf("Solving...\n");
    }
    auto start = chrono::high_resolution_clock::now();
    for (int k = 0; k <= Nt; k += 1){
        this->computeNext(k);
        this->applyReynoldsBC(u_current, u_next);
        this->applyAbsorptionBC();
        // this->computeNextInverted(k + 1);
        // this->applyReynoldsBC(u_next, u_current);
        // this->applyAbsorptionBC();
        if (k % frameRate == 0){
            string fileName = this->outputFolder + "output_data_" + to_string(rank) + "_" + to_string(k/frameRate) + ".vti";
            // saveVTI(u_current, fileName, "Amplitude");
            writeVTI(u_current, fileName, "Amplitude");

        }
        VecSwap(this->u_current, this->u_next);
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

void Solver2d::readInputVtkImageData(string input_file){
    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    printf("Rank: %d\n", rank);
    int* extent = new int[6];
    double* spacing = new double[3];
    vtkSmartPointer<vtkDataArray> data;

    vtkSmartPointer<vtkXMLImageDataReader> reader = vtkSmartPointer<vtkXMLImageDataReader>::New();
    reader->SetFileName(input_file.c_str());
    reader->Update();
    vtkSmartPointer<vtkImageData> imageData = reader->GetOutput();
    extent = imageData->GetExtent();
    spacing = imageData->GetSpacing();
    data = imageData->GetPointData()->GetScalars();

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
    VecSet(this->vel, 3000.0);
    VecAssemblyBegin(this->vel);
    VecAssemblyEnd(this->vel);

    PetscInt NxStart, NzStart, NxEnd, NzEnd;
    DMDAGetCorners(this->da, &NxStart, &NzStart, NULL, &NxEnd, &NzEnd, NULL);

    DMGlobalToLocalBegin(da, this->vel, INSERT_VALUES, this->vel_local);
    DMGlobalToLocalEnd(da, this->vel, INSERT_VALUES, this->vel_local);

    PetscScalar **velArray;
    DMDAVecGetArray(da, this->vel_local, &velArray);

    // PetscInt NzStart, NxStart, NzEnd, NxEnd;
    // DMDAGetCorners(da, &NxStart, &NzStart, NULL, &NxEnd, &NzEnd, NULL);
    for (int j = NzStart; j < NzStart + NzEnd; j++){
        for (int i = NxStart; i < NxStart + NxEnd; i++){
            velArray[j][i] = (PetscScalar)data->GetTuple1(j*this->Nx + i);
        }
    }
    DMDAVecRestoreArray(da, this->vel_local, &velArray);

    DMLocalToGlobalBegin(da,this->vel_local,INSERT_VALUES,this->vel);
    DMLocalToGlobalEnd(da,this->vel_local,INSERT_VALUES,this->vel);
}

void Solver2d::saveVTI(Vec grid, string outputPath, string info){
    PetscInt NxStart, NzStart, NxEnd, NzEnd;
    DMDAGetCorners(da, &NxStart, &NzStart, NULL, &NxEnd, &NzEnd, NULL);
    vtkSmartPointer<vtkImageData> imageData = vtkSmartPointer<vtkImageData>::New();
    imageData->SetExtent(NxStart, NxStart + NxEnd, NzStart, NzStart + NzEnd, 0, 0);
    imageData->SetOrigin(0, 0, 0);
    imageData->SetSpacing(dx, dz, 0);
    vtkSmartPointer<vtkFloatArray> data = vtkSmartPointer<vtkFloatArray>::New();
    data->SetNumberOfComponents(1);
    data->SetNumberOfTuples(NxEnd*NzEnd);
    data->SetName(info.c_str());
    PetscScalar **array;
    DMDAVecGetArray(da, grid, &array);
    // for (int i = NxStart; i < NxStart + NxEnd; i++){
    //     for (int j = NzStart; j < NzStart + NzEnd; j++){
    //         data->SetTuple1(i*this->Nx + j, array[i*this->Nx + j]);
    //     }
    // }
    int jj = 0;
    for (int j = NzStart; j < NzStart + NzEnd; j++){
        int ii = 0;
        for (int i = NxStart; i < NxStart + NxEnd; i++){
            data->SetTuple1(jj*NzEnd + ii, array[j][i]);
            ii++;
        }
        jj++;
    }
    DMDAVecRestoreArray(da, grid, &array);
    imageData->GetPointData()->SetScalars(data);
    vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();
    writer->SetDataModeToAscii();
    writer->SetFileName(outputPath.c_str());
    writer->SetInputData(imageData);
    writer->Write();
}

void Solver2d::writeVTI(Vec grid, string outputPath, string info){
    // printf("Writing VTI file...\n");
    PetscInt NxStart, NzStart, NxEnd, NzEnd;
    DMDAGetCorners(this->da, &NxStart, &NzStart, NULL, &NxEnd, &NzEnd, NULL);
    fstream file;
    file.open(outputPath, ios::out);
    file << "<?xml version=\"1.0\"?>\n";
    file << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    file << "<ImageData WholeExtent=\"" << NxStart << " " << NxStart + NxEnd - 1 << " " << NzStart << " " << NzStart + NzEnd - 1 << " 0 0\" Origin=\"0 0 0\" Spacing=\"" << dx << " " << dz << " 0\">\n";
    file << "<Piece Extent=\"" << NxStart << " " << NxStart + NxEnd - 1 << " " << NzStart << " " << NzStart + NzEnd - 1 << " 0 0\">\n";
    file << "<PointData Scalars=\"" << info << "\">\n";
    file << "<DataArray type=\"Float32\" Name=\"" << info << "\" format=\"ascii\">\n";
    PetscScalar **array;
    DMDAVecGetArray(da, grid, &array);
    for (int j = NzStart; j < NzStart + NzEnd; j++){
        for (int i = NxStart; i < NxStart + NxEnd; i++){
            file << array[j][i] << " ";
        }
        file << "\n";
    }
    DMDAVecRestoreArray(da, grid, &array);
    file << "</DataArray>\n";
    file << "</PointData>\n";
    file << "</Piece>\n";
    file << "</ImageData>\n";
    file << "</VTKFile>\n";
    file.close();
}

void Solver2d::writePVTI(string outputPath, string info){
    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    PetscInt NxStart, NzStart, NxEnd, NzEnd;
    DMDAGetCorners(this->da, &NxStart, &NzStart, NULL, &NxEnd, &NzEnd, NULL);
    if (rank == 0){
        // write pvti file
        string pvtiPath = outputPath.substr(0, outputPath.size() - 4) + ".pvti";
        fstream pvtiFile;
        pvtiFile.open(pvtiPath, ios::out);
        pvtiFile << "<?xml version=\"1.0\"?>\n";
        pvtiFile << "<VTKFile type=\"PImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
        pvtiFile << "<PImageData WholeExtent=\"" << NxStart << " " << NxStart + NxEnd - 1 << " " << NzStart << " " << NzStart + NzEnd - 1 << " 0 0\" Origin=\"0 0 0\" Spacing=\"" << dx << " " << dz << " 0\">\n";
        pvtiFile << "<PPointData Scalars=\"" << info << "\">\n";
        pvtiFile << "<PDataArray type=\"Float32\" Name=\"" << info << "\" format=\"ascii\"/>\n";
        pvtiFile << "</PPointData>\n";
        pvtiFile << "<Piece Extent=\"" << NxStart << " " << NxStart + NxEnd - 1 << " " << NzStart << " " << NzStart + NzEnd - 1 << " 0 0\" Source=\"" << outputPath << "\"/>\n";
        pvtiFile << "</PImageData>\n";
        pvtiFile << "</VTKFile>\n";
        pvtiFile.close();
    }
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
    DMGlobalToLocalBegin(da, this->vel, INSERT_VALUES, this->vel_local);
    DMGlobalToLocalEnd(da, this->vel, INSERT_VALUES, this->vel_local);

    PetscScalar **u_currentArray, **u_nextArray, **velArray;
    DMDAVecGetArray(da, u_current_local, &u_currentArray);
    DMDAVecGetArray(da, u_next_local, &u_nextArray);
    DMDAVecGetArray(da, this->vel_local, &velArray);
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
    bool bc[4];
    NzEnd = NzStart + NzEnd;
    NxEnd = NxStart + NxEnd;
    bc[0] = NxStart == 0;      // left
    bc[1] = NzStart == 0;      // bottom
    bc[2] = NxEnd == this->Nx; // right
    bc[3] = NzEnd == this->Nz; // top
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
    if (bc[0]){
        for (int j = NzStart; j < NzEnd; j++){
            // float courantNumber = dt * vel->get(j, STENCIL)/dx;
            // (*u_next)(j, STENCIL - 1) = (*u_current)(j, STENCIL - 1) - courantNumber*((*u_current)(j, STENCIL) - (*u_current)(j, STENCIL - 1));
            int i = STENCIL + 1;
            float courantNumber = dt * velArray[j][i]/dx;
            u_nextArray[j][i] = u_currentArray[j][i] + courantNumber*(u_currentArray[j][i + 1] - u_currentArray[j][i]);
        }
    }
    if (bc[1]){
        for (int i = NxStart; i < NxEnd; i++){
            // float courantNumber = dt * vel->get(STENCIL, i)/dx;
            // (*u_next)(STENCIL - 1, i) = (*u_current)(STENCIL - 1, i) - courantNumber*((*u_current)(STENCIL, i) - (*u_current)(STENCIL - 1, i));
            int j = STENCIL + 1;
            float courantNumber = dt * velArray[j][i]/dx;
            u_nextArray[j][i] = u_currentArray[j][i] + courantNumber*(u_currentArray[j + 1][i] - u_currentArray[j][i]);
        }
    }
    if (bc[2]){
        for (int j = NzStart; j < NzEnd; j++){
            // float courantNumber = dt * vel->get(j, Nx - STENCIL - 1)/dx;
            // (*u_next)(j, Nx - STENCIL) = (*u_current)(j, Nx - STENCIL) - courantNumber*((*u_current)(j, Nx - STENCIL) - (*u_current)(j, Nx - STENCIL - 1));
            int i = Nx - STENCIL - 2;
            float courantNumber = dt * velArray[j][i]/dx;
            u_nextArray[j][i] = u_currentArray[j][i] - courantNumber*(u_currentArray[j][i] - u_currentArray[j][i - 1]);
        }
    }
    if (bc[3]){
        for (int i = NxStart; i < NxEnd; i++){
            // float courantNumber = dt * vel->get(Nz - STENCIL - 1, i)/dx;
            // (*u_next)(Nz - STENCIL, i) = (*u_current)(Nz - STENCIL, i) - courantNumber*((*u_current)(Nz - STENCIL, i) - (*u_current)(Nz - STENCIL - 1, i));
            int j = Nz - STENCIL - 2;
            float courantNumber = dt * velArray[j][i]/dx;
            u_nextArray[j][i] = u_currentArray[j][i] - courantNumber*(u_currentArray[j][i] - u_currentArray[j - 1][i]);
        }
    }
    DMDAVecRestoreArray(da, u_current_local, &u_currentArray);
    DMDAVecRestoreArray(da, u_next_local, &u_nextArray);

    DMLocalToGlobalBegin(da,u_current_local,INSERT_VALUES,u_current);
    DMLocalToGlobalEnd(da,u_current_local,INSERT_VALUES,u_current);
    DMLocalToGlobalBegin(da,u_next_local,INSERT_VALUES,u_next);
    DMLocalToGlobalEnd(da,u_next_local,INSERT_VALUES,u_next);
}

float Solver2d::mitigation(float x, int border){
    float fat = 0.0055f;
    return expf(-(powf(fat*(border - x), 2.0f)));
}

void Solver2d::applyAbsorptionBC(){
    int border = 25;

    DMGlobalToLocalBegin(da,this->u_current,INSERT_VALUES,this->u_current_local);
    DMGlobalToLocalEnd(da,this->u_current,INSERT_VALUES,this->u_current_local);
    DMGlobalToLocalBegin(da,this->u_next,INSERT_VALUES,this->u_next_local);
    DMGlobalToLocalEnd(da,this->u_next,INSERT_VALUES,this->u_next_local);

    PetscScalar **u_currentArray, **u_nextArray;
    DMDAVecGetArray(da, u_current_local, &u_currentArray);
    DMDAVecGetArray(da, u_next_local, &u_nextArray);

    PetscInt NzStart, NxStart, NzEnd, NxEnd;
    DMDAGetCorners(da, &NxStart, &NzStart, NULL, &NxEnd, &NzEnd, NULL);
    bool bc[4];
    NzEnd = NzStart + NzEnd;
    NxEnd = NxStart + NxEnd;
    bc[0] = NxStart == 0;      // left
    bc[1] = NzStart == 0;      // bottom
    bc[2] = NxEnd == this->Nx; // right
    bc[3] = NzEnd == this->Nz; // top
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

    if (bc[0]){
        for (int j = NzStart; j < NzEnd; j++){
            for (int i = NxStart; i <= border; i++){
                u_currentArray[j][i] = u_currentArray[j][i]*mitigation(i, border);
                u_nextArray[j][i] = u_nextArray[j][i]*mitigation(i, border);
            }
        }
    }
    if (bc[1]){
        for (int i = NxStart; i < NxEnd; i++){
            for (int j = NzStart; j <= border; j++){
                u_currentArray[j][i] = u_currentArray[j][i]*mitigation(j, border);
                u_nextArray[j][i] = u_nextArray[j][i]*mitigation(j, border);
            }
        }
    }
    if (bc[2]){
        for (int j = NzStart; j < NzEnd; j++){
            for (int i = NxEnd - border; i < NxEnd; i++){
                u_currentArray[j][i] = u_currentArray[j][i]*mitigation(Nx - i, border);
                u_nextArray[j][i] = u_nextArray[j][i]*mitigation(Nx - i, border);
            }
        }
    }
    if (bc[3]){
        for (int i = NxStart; i < NxEnd; i++){
            for (int j = NzEnd - border; j < NzEnd; j++){
                u_currentArray[j][i] = u_currentArray[j][i]*mitigation(Nz - j, border);
                u_nextArray[j][i] = u_nextArray[j][i]*mitigation(Nz - j, border);
            }
        }
    }
    DMDAVecRestoreArray(da, u_current_local, &u_currentArray);
    DMDAVecRestoreArray(da, u_next_local, &u_nextArray);

    DMLocalToGlobalBegin(da,u_current_local,INSERT_VALUES,u_current);
    DMLocalToGlobalEnd(da,u_current_local,INSERT_VALUES,u_current);
    DMLocalToGlobalBegin(da,u_next_local,INSERT_VALUES,u_next);
    DMLocalToGlobalEnd(da,u_next_local,INSERT_VALUES,u_next);
}
