#include "Solver2d.h"

Solver2d::Solver2d(string nomeDoArquivo){

    // le os parametros do modelo
    leParametros(nomeDoArquivo);

    // incializa os vetores;
    this->u_current = new Grid2d(this->d.Nz, this->d.Nx);
    this->u_next = new Grid2d(this->d.Nz, this->d.Nx);
    this->sis = new Grid2d(this->d.Nt, this->d.Nx);

    this->posReceptor = this->d.zs/this->d.dz;
    this->modk = 100;

}

Solver2d::~Solver2d(){

}

void Solver2d::imprimeParametros(){

    cout << "Nx = " << this->d.Nx << endl;
    cout << "Nz = " << this->d.Nz << endl;
    cout << "Nt = " << this->d.Nt << endl;

}

void Solver2d::leParametros(string nome){

    // * funcao que le os dados do dominio 2d e a matriz de velocidades

    ifstream myfile;

    myfile.open("../" + nome);

    if(myfile.is_open()){

        // dimensoes do dominio
        myfile >> this->d.X;
        myfile >> this->d.Z;
        myfile >> this->d.T;

        // largura da malha
        myfile >> this->d.dx;
        myfile >> this->d.dz;
        myfile >> this->d.dt;

        // posicao da fonte
        myfile >> this->d.fcorte;
        myfile >> this->d.xs;
        // myfile >> this->d.zs;
        this->d.zs = d.Z/2;

        // numero de iteracoes
        this->d.Nx = this->d.X/this->d.dx;
        this->d.Nz = this->d.Z/this->d.dz;
        this->d.Nt = this->d.T/this->d.dt;

        this->d.vel = new Grid2d(this->d.Nz, this->d.Nx);
        int v;

        // matriz de velocidades
        for (int j = 0; j < this->d.Nz; j++){
            for (int i = 0; i < this->d.Nx; i++){
                // myfile >> v;
                // this->d.vel->set(j, i, v);
                this->d.vel->set(j, i, 2200);
            }
        }

        myfile.close();

    } else {
        cerr << "Falha ao abrir arquivo de parametros" << endl;
        exit;
    }
    
}

void Solver2d::salvaVTI(Dominio d, Grid2d* grid, string nomeDoArq, string info){

    // * Função que gera um arquivo vtk ImageData para o ParaView

    ofstream myfile;

    // cout << "Gerando arquivo data" << to_string(k/modk) << ".vti" << "..." << endl;

    myfile.open("../" + nomeDoArq + ".vti");

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
        cout << "Erro na gravação do arquivo " << nomeDoArq << ".vti" << endl;
    }

}

void Solver2d::salvaVTIbin(Dominio d, Grid2d* grid, string nomeDoArq, string info){

    // * Função que gera um arquivo vtk ImageData para o ParaView

    ofstream myfile;

    // cout << "Gerando arquivo data" << to_string(k/modk) << ".vti" << "..." << endl;

    myfile.open("../" + nomeDoArq + ".vti", ios::binary);

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
        cout << "Erro na gravação do arquivo " << nomeDoArq << ".vti" << endl;
    }

}

void salvaGNUPlot(Dominio d, int k, int modk, string base, Grid2d* u){

    // * Função que gera um arquivo de dados no formato aceito pelo GNUPlot

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

void salvaBIN(Grid2d sis, Dominio d){

    // * grava binario do sismograma
    ofstream file;

    file.open("../sis.bin", ios::out | ios::binary);

    for (int i = 0; i < d.Nx; i++){
        for (int k = 0; k < d.Nt; k++){
            file.write((char *) &sis(k, i), sizeof(float));
        }
    }

    file.close();

}

float Solver2d::fonte(int x, int z, float k){

    // *funcao que simula um pulso sismico na posicao (xs, zs)

    if (x*this->d.dx != this->d.xs || z*this->d.dz != this->d.zs || k*this->d.dt > 0.5){
        return 0;
    } 

    float td = k*this->d.dt - ((2*sqrt(M_PI))/this->d.fcorte);  
    float fc = (this->d.fcorte/(3*sqrt(M_PI)));

    return (1.0 - 2.0 * M_PI * pow(M_PI * fc * td, 2))/pow(M_E, M_PI*pow((M_PI*fc*td), 2));

}

void Solver2d::mdf(Dominio d, Grid2d* u_current, Grid2d* u_next, int k){


    // * Função que calcula o u_next pelo u_current

    float val, courantNumber, const1, const2;

    #pragma omp parallel for collapse(1) private(val, courantNumber, const1, const2)
    for (int j = STENCIL; j <= d.Nz - STENCIL; j++){
        
        for (int i = STENCIL; i <= d.Nx - STENCIL; i++){

            courantNumber = d.dt*d.vel->get(j, i)/d.dx;
            const1 = (pow(courantNumber, 2)/12.0);
            const2 = pow(d.vel->get(j, i)*d.dt, 2);

            val = 
            const1 *
            (
                -1*(u_current->get(j, i - 2) + u_current->get(j - 2, i)) + 
                16*(u_current->get(j, i - 1) + u_current->get(j - 1, i)) - 
                60* u_current->get(j,     i) +
                16*(u_current->get(j, i + 1) + u_current->get(j + 1, i)) -
                   (u_current->get(j, i + 2) + u_current->get(j + 2, i)) 
            ) 
            + 2*u_current->get(j, i) - u_next->get(j, i) - const2 * this->fonte(i, j, k);

            u_next->set(j, i, val);
            
        }
    }

}

void Solver2d::aplicaReynolds(Grid2d* u_current, Grid2d* u_next){

    // * Função que aplica a condição de contorno não-reflexiva de Reynolds

    float courantNumber;

    // * borda esquerda
    //   du/dt - vel*du/dx = 0
    // (u(x,t+dt) - u(x,t))/dt = vel*(u(x+dx,t) - u(x,t))/dx
    // (u(x,t+dt) - u(x,t)) =   dt*(vel*(u(x+dx,t) - u(x,t))/dx
    // u(x, t+dt) = u(x,t) + cou * (u(x+dx,t) - u(x,t) )
    //    onde cou = dt*vel/dx
    #pragma omp parallel for private(courantNumber)
    for(int j = 0; j < d.Nz; j++) {

        courantNumber = d.dt * d.vel->get(j, STENCIL)/d.dx;

        u_next->set(j, STENCIL, u_current->get(j, STENCIL) + courantNumber*(u_current->get(j,STENCIL+1) - u_current->get(j, STENCIL)));
    
    }

    // * borda direita
    //   du/dt + vel*du/dx = 0
    // u(x, t+dt) = u(x,t) - cou * (u(x,t) - u(x-dt,t) )
    // Aqui usamos diferencas atrasadas para discretizar do espaco
    #pragma omp parallel for private(courantNumber)
    for(int j = 0; j < d.Nz; j++) {

        courantNumber = d.dt * d.vel->get(j, d.Nx - STENCIL)/d.dx;

        u_next->set(j, d.Nx-STENCIL, u_current->get(j, d.Nx-STENCIL) - courantNumber*(u_current->get(j, d.Nx-STENCIL) - u_current->get(j, d.Nx-STENCIL-1)));
    
    }

    // * borda superior
    //   du/dt - vel*du/dz = 0
    // u(z, t+dt) = u(z,t) + cou * (u(z+dz,t) - u(z,t) )
    #pragma omp parallel for private(courantNumber)
    for (int i = 0; i < d.Nx; i++) {

        courantNumber = d.dt * d.vel->get(STENCIL, i)/d.dx;

        // u_next->set(STENCIL, i, u_current->get(STENCIL, i) + courantNumber*(u_current->get(STENCIL+1, i) - u_current->get(STENCIL, i)));
        u_next->set(STENCIL, i, u_current->get(STENCIL, i) + courantNumber*(u_current->get(STENCIL+1, i) - u_current->get(STENCIL, i)));
    
    }

    // * borda inferior
    //   du/dt + vel*du/dz = 0
    // u(z, t+dt) = u(z,t) - cou * (u(z,t) - u(z-dz,t) )
    #pragma omp parallel for private(courantNumber)
    for (int i = 0; i < d.Nx; i++) {

        courantNumber = d.dt * d.vel->get(d.Nz - STENCIL, i)/d.dx;

        u_next->set(d.Nz - STENCIL, i, u_current->get(d.Nz-STENCIL, i) - courantNumber*(u_current->get(d.Nz-STENCIL, i) - u_current->get(d.Nz-STENCIL - 1, i)));
    
    }

}

float Solver2d::atenuacao(float x, int borda){

    // * Função aplicada nas bordas para reduzir a amplitude da onda

    float fat = 0.0055;
    return exp(-(pow(fat*(borda - x), 2)));

}

void Solver2d::aplicaAmortecimento(){

    // * Percorre as bordas de das matrizes atual e proxima aplicando uma função de atenuação

    int borda = 23;

    // percorre a faixa superior
    #pragma omp parallel for collapse(1)
    for(int j = STENCIL; j < borda; j++){
        for(int i = STENCIL; i <= d.Nx - STENCIL; i++){

            u_current->set(j, i, u_current->get(j, i)*atenuacao(j, borda));
            u_next->set(j, i, u_next->get(j, i)*atenuacao(j, borda));
            
        }
    }

    // percorre a faixa inferior
    #pragma omp parallel for collapse(1)
    for(int j = d.Nz - borda; j <= d.Nz - STENCIL; j++){
        for(int i = STENCIL; i <= d.Nx - STENCIL; i++){

            u_current->set(j, i, u_current->get(j, i)*atenuacao(d.Nz - j, borda));
            u_next->set(j, i, u_next->get(j, i)*atenuacao(d.Nz - j, borda));
            
        }
    }

    // percorre a faixa esquerda
    #pragma omp parallel for collapse(1)
    for(int j = borda; j < d.Nz - borda; j++){
        for(int i = STENCIL; i <= borda; i++){

            u_current->set(j,i, u_current->get(j, i)*atenuacao(i, borda));
            u_next->set(j,i, u_next->get(j, i)*atenuacao(i, borda));
            
        }
    }

    // percorre a faixa direita
    #pragma omp parallel for collapse(1)
    for(int j = borda; j < d.Nz - borda; j++){
        for(int i = d.Nx - borda; i <= d.Nx - STENCIL; i++){

            u_current->set(j,i, u_current->get(j, i)*atenuacao((d.Nx - i), borda));
            u_next->set(j,i, u_next->get(j, i)*atenuacao((d.Nx - i), borda));
            
        }
    }


}

void Solver2d::solve(){

    // gera arquivo para visualização do modelo de velocidades
    salvaVTI(this->d, this->d.vel, "modelo_velocidades", "Velocidade");

    cout << "Nx = " << d.Nx << endl;
    cout << "Nz = " << d.Nz << endl;
    cout << "Nt = " << d.Nt << endl;

    auto inicio = chrono::high_resolution_clock::now();

    for (int k = 0; k <= d.Nt; k += 2){

        // * calcula u_next
        this->mdf(d, u_current, u_next, k);
        this->aplicaReynolds(u_current, u_next);
        this->aplicaAmortecimento();

        // * armazena na matriz do sismograma
        for (int i = 0; i < d.Nx; i++){
            sis->set(k, i, u_next->get(this->posReceptor, i));
        }
        
        // * calcula u_current
        this->mdf(d, u_next, u_current, k + 1);
        this->aplicaReynolds(u_next, u_current);
        this->aplicaAmortecimento();

        for (int i = 0; i < d.Nx; i++){
            sis->set(k + 1, i, u_current->get(this->posReceptor, i));
        }

        // * gera arquivo de dados a cada 100 iteracoes em k
        if (k % modk == 0){
            string nomeDoArq = "data" + to_string(k/modk);
            salvaVTI(d, u_current, nomeDoArq, "Amplitude");
        }

    }

    salvaBIN(*sis, d);

    cout << "\nArquivos gerados com sucesso!" << endl;

    auto final = chrono::high_resolution_clock::now();
    chrono::duration<double> intervalo = final - inicio;
    cout << "\nTempo decorrido: " << intervalo.count() << "s\n";

}