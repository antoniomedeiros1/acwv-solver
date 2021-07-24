#define _USE_MATH_DEFINES
#define STENCIL 3 

#include <iostream>
#include <fstream>
#include <cmath>
#include "math.h"
#include <string>
#include <chrono>

#include "Grid2d.h"
#include "Dominio.h"

using namespace std;

void leParametros(Dominio* d){

    // * Função que lê os parâmetros da simulação (velocidade constante)

    ifstream myfile;
    myfile.open("../parametros.txt");

    if(myfile.is_open()){

        // dimensoes do dominio
        myfile >> d->X;
        myfile >> d->Z;
        myfile >> d->T;

        // info da onda
        myfile >> d->c;
        myfile >> d->fcorte;

        // largura da malha
        myfile >> d->dx;
        d->dz = d->dx;
        myfile >> d->dt;

        // posicao da fonte
        myfile >> d->xs;
        myfile >> d->zs;

        // constantes
        d->cou = d->dt*d->c/d->dx;
        d->c1 = (pow(d->cou, 2)/12.0); 
        d->c2 = pow(d->c*d->dt, 2);

        // numero de iteracoes
        d->Nx = d->X/d->dx;
        d->Nz = d->Z/d->dz;
        d->Nt = d->T/d->dt;

        myfile.close();

    } else {
        cerr << "Falha ao abrir arquivo de parametros" << endl;
        exit;
    }

}

void leCamposDeVelocidades(Dominio* d){

    // * Função que lê os parâmetros da simulação de um arquivo de velocidades

    ifstream myfile;

    myfile.open("../reservatorio.txt");

    if(myfile.is_open()){

        // dimensoes do dominio
        myfile >> d->X;
        myfile >> d->Z;
        myfile >> d->T;

        // largura da malha
        myfile >> d->dx;
        myfile >> d->dz;
        myfile >> d->dt;

        // posicao da fonte
        myfile >> d->fcorte;
        myfile >> d->xs;
        myfile >> d->zs;

        // numero de iteracoes
        d->Nx = d->X/d->dx;
        d->Nz = d->Z/d->dz;
        d->Nt = d->T/d->dt;

        d->vel = new Grid2d(d->Nz, d->Nx);
        int v;

        // matriz de velocidades
        for (int j = 0; j < d->Nz; j++){
            for (int i = 0; i < d->Nx; i++){
                myfile >> v;
                d->vel->set(j, i, v);
            }
        }

        myfile.close();

    } else {
        cerr << "Falha ao abrir arquivo de parametros" << endl;
        exit;
    }

}

float fonte(int x, int z, float k, Dominio d){

    // *funcao que simula um pulso sismico na posicao (xs, zs)

    if (x*d.dx != d.xs || z*d.dz != d.zs || k*d.dt > 0.5){
        return 0;
    } 

    float td = k*d.dt - ((2*sqrt(M_PI))/d.fcorte);  
    float fc = (d.fcorte/(3*sqrt(M_PI)));

    return (1.0 - 2.0 * M_PI * pow(M_PI * fc * td, 2))/pow(M_E, M_PI*pow((M_PI*fc*td), 2));

}

void mdf(Dominio d, Grid2d* u_current, Grid2d* u_next, int k){

    // * Função que calcula o u_next pelo u_current

    float val, courantNumber, const1, const2;

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
            + 2*u_current->get(j, i) - u_next->get(j, i) - const2 * fonte(i, j, k, d);

            u_next->set(j, i, val);
            
        }
    }

}

void Reynolds(Dominio d, Grid2d* u_current, Grid2d* u_next){

    // * Função que aplica a condição de contorno não-reflexiva de 

    float courantNumber;

    // * borda esquerda
    //   du/dt - vel*du/dx = 0
    // (u(x,t+dt) - u(x,t))/dt = vel*(u(x+dx,t) - u(x,t))/dx
    // (u(x,t+dt) - u(x,t)) =   dt*(vel*(u(x+dx,t) - u(x,t))/dx
    // u(x, t+dt) = u(x,t) + cou * (u(x+dx,t) - u(x,t) )
    //    onde cou = dt*vel/dx
    for(int j = 0; j < d.Nz; j++) {

        courantNumber = d.dt * d.vel->get(j, STENCIL)/d.dx;

        u_next->set(j, STENCIL, u_current->get(j, STENCIL) + courantNumber*(u_current->get(j,STENCIL+1) - u_current->get(j, STENCIL)));
    
    }

    // * borda direita
    //   du/dt + vel*du/dx = 0
    // u(x, t+dt) = u(x,t) - cou * (u(x,t) - u(x-dt,t) )
    // Aqui usamos diferencas atrasadas para discretizar do espaco
    for(int j = 0; j < d.Nz; j++) {

        courantNumber = d.dt * d.vel->get(j, d.Nx - STENCIL)/d.dx;

        u_next->set(j, d.Nx-STENCIL, u_current->get(j, d.Nx-STENCIL) - courantNumber*(u_current->get(j, d.Nx-STENCIL) - u_current->get(j, d.Nx-STENCIL-1)));
    
    }

    // * borda superior
    //   du/dt - vel*du/dz = 0
    // u(z, t+dt) = u(z,t) + cou * (u(z+dz,t) - u(z,t) )
    for (int i = 0; i < d.Nx; i++) {

        courantNumber = d.dt * d.vel->get(STENCIL, i)/d.dx;

        u_next->set(STENCIL, i, u_current->get(STENCIL, i) + courantNumber*(u_current->get(STENCIL+1, i) - u_current->get(STENCIL, i)));
    
    }

    // * borda inferior
    //   du/dt + vel*du/dz = 0
    // u(z, t+dt) = u(z,t) - cou * (u(z,t) - u(z-dz,t) )
    for (int i = 0; i < d.Nx; i++) {

        courantNumber = d.dt * d.vel->get(d.Nz - STENCIL, i)/d.dx;

        u_next->set(d.Nz - STENCIL, i, u_current->get(d.Nz-STENCIL, i) - courantNumber*(u_current->get(d.Nz-STENCIL, i) - u_current->get(d.Nz-STENCIL - 1, i)));
    
    }

}

float atenuacao(float x, int borda){

    // * Função aplicada nas bordas para reduzir a amplitude da onda

    float fat = 0.0035;
    return exp(-(pow(fat*(borda - x), 2)));

}

void camadasDeAbsorcao(Dominio d, Grid2d* u_current, Grid2d* u_next){

    // * Percorre as bordas de das matrizes atual e proxima aplicando uma função de atenuação

    int borda = 23;

    // percorre a faixa superior
    for(int j = STENCIL; j < borda; j++){
        for(int i = STENCIL; i <= d.Nx - STENCIL; i++){

            u_current->set(j, i, u_current->get(j, i)*atenuacao(j, borda));
            u_next->set(j, i, u_next->get(j, i)*atenuacao(j, borda));
            
        }
    }

    // percorre a faixa inferior
    for(int j = d.Nz - borda; j <= d.Nz - STENCIL; j++){
        for(int i = STENCIL; i <= d.Nx - STENCIL; i++){

            u_current->set(j, i, u_current->get(j, i)*atenuacao(d.Nz - j, borda));
            u_next->set(j, i, u_next->get(j, i)*atenuacao(d.Nz - j, borda));
            
        }
    }

    // percorre a faixa esquerda
    for(int j = borda; j < d.Nz - borda; j++){
        for(int i = STENCIL; i <= borda; i++){

            u_current->set(j,i, u_current->get(j, i)*atenuacao(i, borda));
            u_next->set(j,i, u_next->get(j, i)*atenuacao(i, borda));
            
        }
    }

    // percorre a faixa direita
    for(int j = borda; j < d.Nz - borda; j++){
        for(int i = d.Nx - borda; i <= d.Nx - STENCIL; i++){

            u_current->set(j,i, u_current->get(j, i)*atenuacao((d.Nx - i), borda));
            u_next->set(j,i, u_next->get(j, i)*atenuacao((d.Nx - i), borda));
            
        }
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

void salvaVTI(Dominio d, Grid2d* u, string nomeDoArq, string info){

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
                myfile << u->get(j, i) << " ";
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

int main() {

    // * inicializa os parametros da simulação
    Dominio d;
    leCamposDeVelocidades(&d);

    // * salva a matriz de velocidades em vti
    salvaVTI(d, d.vel, "modelo_vel", "velocidade");

    // * é necessária a utilização de apenas duas matrizes para implementar o método de diferenças finitas
    Grid2d u_current(d.Nz, d.Nx);
    Grid2d u_next(d.Nz, d.Nx);

    // * matriz para o sismograma com uma dimensão no espaço e uma no tempo
    Grid2d sis(d.Nt, d.Nx);
    int posicao_receptor = 30; // profundidade em pontos dos receptores
    
    int modk = 50;

    cout << "Nx = " << d.Nx << endl;
    cout << "Nz = " << d.Nz << endl;
    cout << "Nt = " << d.Nt << endl;

    auto inicio = chrono::high_resolution_clock::now();

    for (int k = 0; k <= d.Nt; k += 2){

        // * calcula u_next
        mdf(d, &u_current, &u_next, k);
        Reynolds(d, &u_current, &u_next);
        camadasDeAbsorcao(d, &u_current, &u_next);

        // * armazena na matriz do sismograma
        for (int i = 0; i < d.Nx; i++){
            sis(k, i) = u_next(posicao_receptor, i);
        }
        
        // * calcula u_current
        mdf(d, &u_next, &u_current, k + 1);
        Reynolds(d, &u_next, &u_current);
        camadasDeAbsorcao(d, &u_next, &u_current);

        for (int i = 0; i < d.Nx; i++){
            sis(k + 1, i) = u_current(posicao_receptor, i);
        }

        // * gera arquivo de dados a cada 100 iteracoes em k
        if (k % modk == 0){
            string nomeDoArq = "data" + to_string(k/modk);
            salvaVTI(d, &u_current, nomeDoArq, "Amplitude");
        }

    }

    // * grava binario do sismograma
    ofstream file;

    file.open("../sis.bin", ios::out | ios::binary);

    for (int i = 0; i < d.Nx; i++){
        for (int k = 0; k < d.Nt; k++){
            file.write((char *) &sis(k, i), sizeof(float));
        }
    }

    file.close();

    cout << "\nArquivos gerados com sucesso!" << endl;

    auto final = chrono::high_resolution_clock::now();
    chrono::duration<double> intervalo = final - inicio;
    cout << "\nTempo decorrido: " << intervalo.count() << "s\n";

    return 0;
}
