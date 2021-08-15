#include "Solver3d.h"

Solver3d::Solver3d(){

    // lendo os parametros
    cout << "Lendo os parametros..." << endl;
    this->leParametros("parametros3d.txt");

    cout << "Parametros lidos com sucesso!" << endl;
    cout << endl;
    cout << "X  = " << this->d.X  << "m" << endl;
    cout << "Y  = " << this->d.Y  << "m" << endl;
    cout << "Z  = " << this->d.Z  << "m" << endl;
    cout << "T  = " << this->d.T  << "s" << endl;
    cout << "dx = " << this->d.dx << "m" << endl;
    cout << "dt = " << this->d.dt << "s" << endl;
    cout << "Nx = " << this->d.Nx << "m" << endl;
    cout << "Ny = " << this->d.Ny << "m" << endl;
    cout << "Nz = " << this->d.Nz << "m" << endl;
    cout << "xs = " << this->d.xs << "m" << endl;
    cout << "ys = " << this->d.ys << "m" << endl;
    cout << "zs = " << this->d.zs << "m" << endl;
    cout << endl;

    cout << "Inicializando vetores..." << endl << endl;
    // inicializa vetores
    this->u_current = new Grid3d(this->d.Nx, this->d.Ny, this->d.Nz);
    this->u_next    = new Grid3d(this->d.Nx, this->d.Ny, this->d.Nz);

}

Solver3d::~Solver3d(){}

void Solver3d::leParametros(string nome){

    // * funcao que le os parâmetros do dominio 3D

    ifstream myfile;

    myfile.open("../" + nome);

    if(myfile.is_open()){

        // dimensoes do dominio
        myfile >> this->d.X;
        myfile >> this->d.Y;
        myfile >> this->d.Z;
        myfile >> this->d.T;

        // largura da malha
        myfile >> this->d.dx;
        this->d.dy = this->d.dx;
        this->d.dz = this->d.dx;
        myfile >> this->d.dt;

        // posicao da fonte
        myfile >> this->d.fcorte;
        myfile >> this->d.xs;
        myfile >> this->d.ys;
        myfile >> this->d.zs;

        // numero de iteracoes
        this->d.Nx = this->d.X/this->d.dx;
        this->d.Ny = this->d.Y/this->d.dy;
        this->d.Nz = this->d.Z/this->d.dz;
        this->d.Nt = this->d.T/this->d.dt;

        // velocidade da onda
        this->d.c = 2200;

        // this->d.vel = new Grid2d(this->d.Nz, this->d.Nx);
        // int v;

        // // matriz de velocidades
        // for (int j = 0; j < this->d.Nz; j++){
        //     for (int i = 0; i < this->d.Nx; i++){
        //         myfile >> v;
        //         this->d.vel->set(j, i, v);
        //     }
        // }

        myfile.close();

    } else {
        cerr << "Falha ao abrir arquivo de parametros" << endl;
        exit;
    }
    
}

void Solver3d::salvaVTI(Dominio d, Grid3d* u, string nomeDoArq, string info){


    // * Função que gera um arquivo vtk ImageData para o ParaView

    ofstream myfile;

    // cout << "Gerando arquivo data" << to_string(k/modk) << ".vti" << "..." << endl;

    myfile.open("../" + nomeDoArq + ".vti");

    if(myfile.is_open()){

        myfile << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
        myfile << "  <ImageData WholeExtent= \"" <<  STENCIL << " " << d.Nx - 1 - STENCIL << " " << STENCIL << " " << d.Ny - 1 - STENCIL << " " << STENCIL << " " << d.Nz - 1 - STENCIL << "\" ";
        myfile << "Origin = \"" << STENCIL << " " << STENCIL << " " << d.Nz - 1 - STENCIL << "\" ";
        myfile << "Spacing = \"" << d.dx << " " << d.dy << " " << d.dz << "\">\n";
        myfile << "    <Piece Extent = \"" << STENCIL << " " << d.Nx - 1 - STENCIL << " " << STENCIL << " " << d.Ny - 1 - STENCIL << " " << STENCIL << " " << d.Nz - 1 - STENCIL << "\">\n";
        myfile << "      <PointData Scalars=\"" + info + "\">\n";
        myfile << "        <DataArray type=\"Float32\" Name=\"" + info + "\" format=\"ascii\">\n";

        for (int k = STENCIL; k < d.Nz - STENCIL; k++){
            for (int j = STENCIL; j < d.Ny - STENCIL; j++){
                for (int i = STENCIL; i < d.Nx - STENCIL; i++){
                    myfile << u->get(k, j, i) << " ";
                }
            }
        }

        myfile << "\n        </DataArray>";
        myfile << "\n      </PointData>";
        myfile << "\n    </Piece>";
        myfile << "\n  </ImageData>";
        myfile << "\n</VTKFile>";

    } else {
        cout << "Erro na gravação do arquivo " << nomeDoArq << ".vti" << endl;
    }


}

// void Solver3d::salvaVTIBinary(Dominio d, Grid3d* u, string nomeDoArq, string info){


//     // * Função que gera um arquivo vtk ImageData para o ParaView

//     ofstream myfile;

//     myfile.open("../" + nomeDoArq + ".vti", ios::binary);

//     if(myfile.is_open()){

//         myfile << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
//         myfile << "  <ImageData WholeExtent= \"" <<  STENCIL << " " << d.Nx - 1 - STENCIL << " " << STENCIL << " " << d.Ny - 1 - STENCIL << " " << STENCIL << " " << d.Nz - 1 - STENCIL << "\" ";
//         myfile << "Origin = \"" << STENCIL << " " << STENCIL << " " << d.Nz - 1 - STENCIL << "\" ";
//         myfile << "Spacing = \"" << d.dx << " " << d.dy << " " << d.dz << "\">\n";
//         myfile << "    <Piece Extent = \"" << STENCIL << " " << d.Nx - 1 - STENCIL << " " << STENCIL << " " << d.Ny - 1 - STENCIL << " " << STENCIL << " " << d.Nz - 1 - STENCIL << "\">\n";
//         myfile << "      <PointData Scalars=\"" + info + "\">\n";
//         myfile << "        <DataArray type=\"Float32\" Name=\"" + info + "\" format=\"binary\">\n";

//         int result_size;
//         char *encoding = base64((char *)u->firstptr(), sizeof(float) * u->getSize(), &result_size);
//         myfile.write(encoding, result_size);

//         myfile << "\n        </DataArray>";
//         myfile << "\n      </PointData>";
//         myfile << "\n    </Piece>";
//         myfile << "\n  </ImageData>";
//         myfile << "\n</VTKFile>";

//     } else {
//         cout << "Erro na gravação do arquivo " << nomeDoArq << ".vti" << endl;
//     }


// }


float Solver3d::pulso3d(int x, int y, int z, int t){

    // *funcao que simula um pulso sismico na posicao (xs, zs)

    if (x*this->d.dx != this->d.xs || y*this->d.dy != this->d.ys || z*this->d.dz != this->d.zs){
        return 0;
    } 

    float td = t*this->d.dt - ((2*sqrt(M_PI))/this->d.fcorte);  
    float fc = (this->d.fcorte/(3*sqrt(M_PI)));
    float r = (1.0 - 2.0 * M_PI * pow(M_PI * fc * td, 2))/pow(M_E, M_PI*pow((M_PI*fc*td), 2));
    // cout << r << endl;
    return r;

}

void Solver3d::mdf3d(Grid3d* u1, Grid3d* u2, int t){

    // * Calcula u2 a partir de u1 por DF

    float val;
    float courantNumber = d.dt*d.c/d.dx;
    float const1 = (pow(courantNumber, 2)/12.0);
    float const2 = pow(d.c*d.dt, 2);

    for (int k = STENCIL; k < d.Nz - STENCIL; k++){
        for (int j = STENCIL; j < d.Ny - STENCIL; j++){
            for (int i = STENCIL; i < d.Nx - STENCIL; i++){

                val = 
                const1 *
                (
                -1*(u1->get(k,j,i-2) + u1->get(k,j,i+2) + u1->get(k,j-2,i) + u1->get(k,j+2,i) + u1->get(k-2,j,i) + u1->get(k+2,j,i)) +
                16*(u1->get(k,j,i-1) + u1->get(k,j,i+1) + u1->get(k,j-1,i) + u1->get(k,j+1,i) + u1->get(k-1,j,i) + u1->get(k+1,j,i)) -
                90*u1->get(k, j, i)
                ) + 
                2*u1->get(k, j, i) - u2->get(k, j, i) + this->pulso3d(i, j, k, t);

                u2->set(k, j, i, val);

            }
        }
    }

}


void Solver3d::mdf3d(Grid3d* u1, Grid3d* u2){

    // * Calcula u2 a partir de u1 por DF

    float val;
    float courantNumber = d.dt*d.c/d.dx;
    float const1 = (pow(courantNumber, 2)/12.0);

    for (int k = STENCIL; k < d.Nz - STENCIL; k++){
        for (int j = STENCIL; j < d.Ny - STENCIL; j++){
            for (int i = STENCIL; i < d.Nx - STENCIL; i++){

                val = 
                const1 *
                (
                -1*(u1->get(k,j,i-2) + u1->get(k,j,i+2) + u1->get(k,j-2,i) + u1->get(k,j+2,i) + u1->get(k-2,j,i) + u1->get(k+2,j,i)) +
                16*(u1->get(k,j,i-1) + u1->get(k,j,i+1) + u1->get(k,j-1,i) + u1->get(k,j+1,i) + u1->get(k-1,j,i) + u1->get(k+1,j,i)) -
                90*u1->get(k, j, i)
                ) + 
                2*u1->get(k, j, i) - u2->get(k, j, i);

                u2->set(k, j, i, val);

            }
        }
    }

}

void Solver3d::solve(){

    int modk = 100;
    int t; // iterador temporal

    auto inicio = chrono::high_resolution_clock::now();

    // laço temporal com fonte
    cout << "Iniciando laço temporal" << endl;
    for (t = 0; t < d.Nt; t += 2){
        
        // * calcula u_next
        this->mdf3d(u_current, u_next, t);

        // * gera arquivo de dados a cada 100 iteracoes em k
        if (t % modk == 0){
            string nomeDoArq = "data" + to_string(t/modk);
            this->salvaVTI(d, u_current, nomeDoArq, "Amplitude");
        }
        
        // * calcula u_current
        this->mdf3d(u_next, u_current, t + 1);

    }

    cout << "\nArquivos gerados com sucesso!" << endl;

    auto final = chrono::high_resolution_clock::now();
    chrono::duration<double> intervalo = final - inicio;
    cout << "\nTempo decorrido: " << intervalo.count() << "s\n";

}