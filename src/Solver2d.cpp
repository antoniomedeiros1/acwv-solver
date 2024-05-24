#include "../include/Solver2d.h"


Solver2d::Solver2d(string input_file, string output_folder, int number_of_steps, float dt, int number_of_frames){
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
    this->d.fcorte = 40;
    this->d.xs = int(this->d.X/2);
    this->d.zs = int(this->d.Z/2);
    this->u_current = new Grid2d(this->d.Nz, this->d.Nx);
    this->u_next    = new Grid2d(this->d.Nz, this->d.Nx);
    this->f         = new Grid2d(this->d.Nz, this->d.Nx);
}

Solver2d::~Solver2d(){}

void Solver2d::solve(){
    printParameters();
    printf("Saving velocity field...\n");
    saveVTI(this->d, this->d.vel, this->outputFolder + "velocity_field", "Velocity");
    printf("Solving...\n");
    auto start = chrono::high_resolution_clock::now();
    for (int k = 0; k <= d.Nt; k += 2){
        this->computeNext(d, u_current, u_next, k);
        this->applyReynoldsBC(u_current, u_next);
        this->applyAbsorptionBC();
        this->computeNext(d, u_next, u_current, k + 1);
        this->applyReynoldsBC(u_next, u_current);
        this->applyAbsorptionBC();
        if (k % frameRate == 0){
            string fileName = this->outputFolder + "output_data" + to_string(k/frameRate);
            saveVTI(d, u_current, fileName, "Amplitude");
        }
    }
    auto final = chrono::high_resolution_clock::now();
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
    myfile.open(input_file);
    if(myfile.is_open()){
        myfile >> this->d.Nx; 
        myfile >> this->d.Nz;
        myfile >> this->d.Nt;
        myfile >> this->d.dx;
        this->d.dz = this->d.dx;
        myfile >> this->d.dt;
        this->d.vel = new Grid2d(this->d.Nz, this->d.Nx);
        float v;
        for (int i = 0; i < this->d.Nx; i++){
            for (int j = 0; j < this->d.Nz; j++){
                myfile >> v;
                this->d.vel->set(j, i, v);
            }
        }
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
    this->d.vel = new Grid2d(this->d.Nz, this->d.Nx);
    for (int i = 0; i < this->d.Nx; i++){
        for (int j = 0; j < this->d.Nz; j++){
            this->d.vel->set(j, i, data->GetTuple1(i + j*this->d.Nx));
        }
    }
}

void Solver2d::saveVTI(Domain d, Grid2d* grid, string outputPath, string info){
    ofstream myfile;
    myfile.open(outputPath + ".vti");
    if(myfile.is_open()){
        myfile << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
        myfile << "  <ImageData WholeExtent= \"" <<  STENCIL << " " << d.Nx - 1 - STENCIL << " " << STENCIL << " " << d.Nz - 1 - STENCIL << " " << 0 << " " << 0 << "\" ";
        myfile << "Origin = \"" << STENCIL << " " << d.Nz - 1 << " " << 0 << "\" ";
        myfile << "Spacing = \"" << d.dx << " " << d.dz << " " << 0 << "\">\n";
        myfile << "    <Piece Extent = \"" << STENCIL << " " << d.Nx - 1 - STENCIL << " " << STENCIL << " " << d.Nz - 1 - STENCIL << " " << 0 << " " << 0 << "\">\n";
        myfile << "      <PointData Scalars=\"" + info + "\">\n";
        myfile << "        <DataArray type=\"Float32\" Name=\"" + info + "\" format=\"ascii\">\n";
        for (int j = STENCIL; j < d.Nz - STENCIL; j++){
            for (int i = STENCIL; i < d.Nx - STENCIL; i++){
                myfile << grid->get(j, i) << " ";
            }
        }
        myfile << "\n        </DataArray>";
        myfile << "\n      </PointData>";
        myfile << "\n    </Piece>";
        myfile << "\n  </ImageData>";
        myfile << "\n</VTKFile>";
    } else {
        cout << "Failed to create " << outputPath << ".vti" << endl;
    }
}

void Solver2d::saveVTIbin(Domain d, Grid2d* grid, string outputPath, string info){
    ofstream myfile;
    myfile.open(outputPath + ".vti", ios::binary);
    if(myfile.is_open()){
        myfile << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
        myfile << "  <ImageData WholeExtent= \"" <<  STENCIL << " " << d.Nx - 1 - STENCIL << " " << STENCIL << " " << d.Nz - 1 - STENCIL << " " << 0 << " " << 0 << "\" ";
        myfile << "Origin = \"" << STENCIL << " " << d.Nz - 1 << " " << 0 << "\" ";
        myfile << "Spacing = \"" << d.dx << " " << d.dz << " " << 0 << "\">\n";
        myfile << "    <Piece Extent = \"" << STENCIL << " " << d.Nx - 1 - STENCIL << " " << STENCIL << " " << d.Nz - 1 - STENCIL << " " << 0 << " " << 0 << "\">\n";
        myfile << "      <PointData Scalars=\"" + info + "\">\n";
        myfile << "        <DataArray type=\"Float32\" Name=\"" + info + "\" format=\"binary\">\n";
        myfile.write((char *) grid->firstptr(), grid->getSize() * sizeof(float));
        myfile << "\n        </DataArray>";
        myfile << "\n      </PointData>";
        myfile << "\n    </Piece>";
        myfile << "\n  </ImageData>";
        myfile << "\n</VTKFile>";
    } else {
        cout << "Erro na gravação do arquivo " << outputPath << ".vti" << endl;
    }
}

void salvaGNUPlot(Domain d, int k, int modk, string base, Grid2d* u){
    ofstream myfile;
    cout << "Gerando arquivo data" << to_string(k/modk) + base << "..." << endl;
    myfile.open("../data" + to_string(k/modk) + base);
        for (int j = STENCIL; j < d.Nz - STENCIL; j++){
            for (int i = STENCIL; i < d.Nx - STENCIL; i++){
                myfile << j*d.dz << " " << i*d.dx << " " << u->get(j, i) << "\n";
            }
            myfile << "\n\n";
        }
    myfile.close();
}

float Solver2d::source(int x, int z, float k){
    if (x != (int)(this->d.xs/this->d.dx) || z!= (int)(this->d.zs/this->d.dz) || k*this->d.dt > 0.5){
        return 0;
    } 
    float td = k*this->d.dt - ((2.0f*sqrtf(M_PI))/this->d.fcorte);  
    float fc = (this->d.fcorte/(3.0f*sqrtf(M_PI)));
    return (1.0f - 2.0f * M_PI * powf(M_PI * fc * td, 2.0f))/powf(M_E, M_PI*powf((M_PI*fc*td), 2.0f));
}

void Solver2d::computeNext(Domain d, Grid2d* u_current, Grid2d* u_next, int k){
    int sizez = d.Nz - STENCIL;
    int sizex = d.Nx - STENCIL;
    float val, courantNumber, const1, const2;
    for (int j = STENCIL; j <= sizez; j++){
        for (int i = STENCIL; i <= sizex; i++){
            courantNumber = d.dt*d.vel->get(j, i)/d.dx;
            const1 = (powf(courantNumber, 2.0f)/12.0f);
            const2 = powf(d.vel->get(j, i)*d.dt, 2.0f);
            val = 
            const1 *
            (
                -1*(u_current->get(j, i - 2) + u_current->get(j - 2, i)) + 
                16*(u_current->get(j, i - 1) + u_current->get(j - 1, i)) - 
                60* u_current->get(j,     i) +
                16*(u_current->get(j, i + 1) + u_current->get(j + 1, i)) -
                   (u_current->get(j, i + 2) + u_current->get(j + 2, i)) 
            ) 
            + 2*u_current->get(j, i) - u_next->get(j, i) - const2 * this->source(i, j, k);
            u_next->set(j, i, val);
        }
    }
}

void Solver2d::applyReynoldsBC(Grid2d* u_current, Grid2d* u_next){
    for(int j = STENCIL; j < d.Nz - STENCIL; j++){
        for(int i = STENCIL; i <= STENCIL + 1; i++){
            float courantNumber = d.dt * d.vel->get(j, i)/d.dx;
            (*u_next)(j, i) = (*u_current)(j, i) + courantNumber*((*u_current)(j, i + 1) - (*u_current)(j, i));
        }
    }
    for(int j = STENCIL; j < d.Nz - STENCIL; j++){
        for(int i = d.Nx - STENCIL - 1; i <= d.Nx - STENCIL; i++){
            float courantNumber = d.dt * d.vel->get(j, i)/d.dx;
            (*u_next)(j, i) = (*u_current)(j, i) - courantNumber*((*u_current)(j, i) - (*u_current)(j, i - 1));
        }
    }
    for(int i = STENCIL; i < d.Nx - STENCIL; i++){
        for(int j = STENCIL; j <= STENCIL + 1; j++){
            float courantNumber = d.dt * d.vel->get(j, i)/d.dx;
            (*u_next)(j, i) = (*u_current)(j, i) + courantNumber*((*u_current)(j + 1, i) - (*u_current)(j, i));
        }
    }
    for(int i = STENCIL; i < d.Nx - STENCIL; i++){
        for(int j = d.Nz - STENCIL - 1; j <= d.Nz - STENCIL; j++){
            float courantNumber = d.dt * d.vel->get(j, i)/d.dx;
            (*u_next)(j, i) = (*u_current)(j, i) - courantNumber*((*u_current)(j, i) - (*u_current)(j - 1, i));
        }
    }
}

float Solver2d::mitigation(float x, int borda){
    float fat = 0.0055f;
    return expf(-(powf(fat*(borda - x), 2.0f)));
}

void Solver2d::applyAbsorptionBC(){
    int borda = 25;
    for(int j = STENCIL; j < borda; j++){
        for(int i = STENCIL; i <= d.Nx - STENCIL; i++){
            u_current->set(j, i, u_current->get(j, i)*mitigation(j, borda));
            u_next->set(j, i, u_next->get(j, i)*mitigation(j, borda));
        }
    }
    for(int j = d.Nz - borda; j <= d.Nz - STENCIL; j++){
        for(int i = STENCIL; i <= d.Nx - STENCIL; i++){
            u_current->set(j, i, u_current->get(j, i)*mitigation(d.Nz - j, borda));
            u_next->set(j, i, u_next->get(j, i)*mitigation(d.Nz - j, borda));
        }
    }
    for(int j = borda; j < d.Nz - borda; j++){
        for(int i = STENCIL; i <= borda; i++){
            u_current->set(j,i, u_current->get(j, i)*mitigation(i, borda));
            u_next->set(j,i, u_next->get(j, i)*mitigation(i, borda));
        }
    }
    for(int j = borda; j < d.Nz - borda; j++){
        for(int i = d.Nx - borda; i <= d.Nx - STENCIL; i++){
            u_current->set(j,i, u_current->get(j, i)*mitigation((d.Nx - i), borda));
            u_next->set(j,i, u_next->get(j, i)*mitigation((d.Nx - i), borda));
        }
    }
}

