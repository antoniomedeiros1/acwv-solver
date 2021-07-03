#define _USE_MATH_DEFINES
#define STENCIL 3 

#include <iostream>
#include <fstream>
#include <cmath>
#include "math.h"
#include <string>
#include <chrono>

#include "Matriz2d.h"
#include "Dominio.h"

using namespace std;

void leParametros(Dominio* d){

    // * Função que lê os parâmetros da simulação

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

float fonte(int x, int z, float k, Dominio d){

    // *funcao que simula um pulso sismico na posicao (xs, zs)

    if (x*d.dx != d.xs || z*d.dz != d.zs || k*d.dt > 0.5){
        return 0;
    } 

    float td = k*d.dt - ((2*sqrt(M_PI))/d.fcorte);  
    float fc = (d.fcorte/(3*sqrt(M_PI)));

    return (1.0 - 2.0 * M_PI * pow(M_PI * fc * td, 2))/pow(M_E, M_PI*pow((M_PI*fc*td), 2));

}

void mdf(Dominio d, Matriz2d* u_current, Matriz2d* u_next, int k){

    // * Função que calcula o u_next pelo u_current

    float val;

    for (int i = STENCIL; i <= d.Nx - STENCIL; i++){
        
        for (int j = STENCIL; j <= d.Nz - STENCIL; j++){

            val = 
            d.c1 *
            (
                -1*(u_current->get(i - 2, j) + u_current->get(i, j - 2)) + 
                16*(u_current->get(i - 1, j) + u_current->get(i, j - 1)) - 
                60* u_current->get(i,     j) +
                16*(u_current->get(i + 1, j) + u_current->get(i, j + 1)) -
                   (u_current->get(i + 2, j) + u_current->get(i, j + 2)) 
            ) 
            + 2*u_current->get(i, j) - u_next->get(i, j) - d.c2 * fonte(i, j, k, d);

            u_next->set(i,j, val);
            
        }
    }

}

void Reynolds(Dominio d, Matriz2d* u_current, Matriz2d* u_next){

    // * Função que aplica a condição de contorno não-reflexiva de 

    // * borda esquerda
    //   du/dt - vel*du/dx = 0
    // (u(x,t+dt) - u(x,t))/dt = vel*(u(x+dx,t) - u(x,t))/dx
    // (u(x,t+dt) - u(x,t)) =   dt*(vel*(u(x+dx,t) - u(x,t))/dx
    // u(x, t+dt) = u(x,t) + cou * (u(x+dx,t) - u(x,t) )
    //    onde cou = dt*vel/dx
    for(int j = 0; j < d.Nz; j++) {
        u_next->set(STENCIL, j, u_current->get(STENCIL,j) + d.cou*(u_current->get(STENCIL+1,j) - u_current->get(STENCIL,j)));
    }

    // * borda direita
    //   du/dt + vel*du/dx = 0
    // u(x, t+dt) = u(x,t) - cou * (u(x,t) - u(x-dt,t) )
    // Aqui usamos diferencas atrasadas para discretizar do espaco
    for(int j = 0; j < d.Nz; j++) {
        u_next->set(d.Nx-STENCIL, j, u_current->get(d.Nx-STENCIL,j) - d.cou*(u_current->get(d.Nx-STENCIL,j) - u_current->get(d.Nx-STENCIL-1,j)));
    }

    // * borda superior
    //   du/dt - vel*du/dz = 0
    // u(z, t+dt) = u(z,t) + cou * (u(z+dz,t) - u(z,t) )
    for (int i = 0; i < d.Nx; i++) {
        u_next->set(i, STENCIL, u_current->get(i,STENCIL) + d.cou*(u_current->get(i,STENCIL+1) - u_current->get(i,STENCIL)));
    }

    // * borda inferior
    //   du/dt + vel*du/dz = 0
    // u(z, t+dt) = u(z,t) - cou * (u(z,t) - u(z-dz,t) )
    for (int i = 0; i < d.Nx; i++) {
        u_next->set(i, d.Nz - STENCIL, u_current->get(i,d.Nz-STENCIL) - d.cou*(u_current->get(i,d.Nz-STENCIL) - u_current->get(i,d.Nz-STENCIL - 1)));
    }

}

float atenuacao(float x, int borda){

    // * Função aplicada nas bordas para reduzir a amplitude da onda

    float fat = 0.00005;
    return exp(-(fat*pow(borda - x, 2)));

}

int dist(int i, int j, Dominio d){

    // * Função que retorna a menor distancia em pontos entre o ponto (i, j) e o contorno 

    return min(min(i, d.Nx - i), min(j, d.Nz - j));

}

void mdfComAbsorcao(Dominio d, Matriz2d* u_current, Matriz2d* u_next, int k){

    // * Função que calcula o u_next pelo u_current

    float val;
    int borda = 20;


    // percorre a faixa superior
    for(int i = STENCIL; i <= d.Nx - STENCIL; i++){
        for(int j = STENCIL; j < borda; j++){

            val = 
            d.c1 *
            (
                -1*(u_current->get(i - 2, j) + u_current->get(i, j - 2)) + 
                16*(u_current->get(i - 1, j) + u_current->get(i, j - 1)) - 
                60* u_current->get(i,     j) +
                16*(u_current->get(i + 1, j) + u_current->get(i, j + 1)) -
                   (u_current->get(i + 2, j) + u_current->get(i, j + 2)) 
            ) 
            + 2*u_current->get(i, j) - u_next->get(i, j);

            val *= atenuacao(j, borda);

            u_next->set(i,j, val);
            
        }
    }

    // percorre a faixa esquerda
    for(int i = STENCIL; i < borda; i++){
        for(int j = borda; j < d.Nz - borda; j++){

            val = 
            d.c1 *
            (
                -1*(u_current->get(i - 2, j) + u_current->get(i, j - 2)) + 
                16*(u_current->get(i - 1, j) + u_current->get(i, j - 1)) - 
                60* u_current->get(i,     j) +
                16*(u_current->get(i + 1, j) + u_current->get(i, j + 1)) -
                   (u_current->get(i + 2, j) + u_current->get(i, j + 2)) 
            ) 
            + 2*u_current->get(i, j) - u_next->get(i, j);

            val *= atenuacao(i, borda);

            u_next->set(i,j, val);
            
        }
    }

    // percorre o miolo
    for(int i = borda; i < d.Nx - borda; i++){
        for(int j = borda; j < d.Nz - borda; j++){

            val = 
            d.c1 *
            (
                -1*(u_current->get(i - 2, j) + u_current->get(i, j - 2)) + 
                16*(u_current->get(i - 1, j) + u_current->get(i, j - 1)) - 
                60* u_current->get(i,     j) +
                16*(u_current->get(i + 1, j) + u_current->get(i, j + 1)) -
                   (u_current->get(i + 2, j) + u_current->get(i, j + 2)) 
            ) 
            + 2*u_current->get(i, j) - u_next->get(i, j) - d.c2 * fonte(i, j, k, d);

            u_next->set(i,j, val);
            
        }
    }

    // percorre a faixa direita
    for(int i = d.Nx - borda; i <= d.Nx - STENCIL; i++){
        for(int j = borda; j < d.Nz - borda; j++){

            val = 
            d.c1 *
            (
                -1*(u_current->get(i - 2, j) + u_current->get(i, j - 2)) + 
                16*(u_current->get(i - 1, j) + u_current->get(i, j - 1)) - 
                60* u_current->get(i,     j) +
                16*(u_current->get(i + 1, j) + u_current->get(i, j + 1)) -
                   (u_current->get(i + 2, j) + u_current->get(i, j + 2)) 
            ) 
            + 2*u_current->get(i, j) - u_next->get(i, j);

            val *= atenuacao((d.Nx - i), borda);

            u_next->set(i,j, val);
            
        }
    }
    
    // percorre a faixa inferior
    for(int i = STENCIL; i <= d.Nx - STENCIL; i++){
        for(int j = d.Nz - borda; j <= d.Nz - STENCIL; j++){

            val = 
            d.c1 *
            (
                -1*(u_current->get(i - 2, j) + u_current->get(i, j - 2)) + 
                16*(u_current->get(i - 1, j) + u_current->get(i, j - 1)) - 
                60* u_current->get(i,     j) +
                16*(u_current->get(i + 1, j) + u_current->get(i, j + 1)) -
                   (u_current->get(i + 2, j) + u_current->get(i, j + 2)) 
            ) 
            + 2*u_current->get(i, j) - u_next->get(i, j);

            val *= atenuacao(d.Nz - j, borda);

            u_next->set(i,j, val);
            
        }
    }

}

void salvaGNUPlot(Dominio d, int k, int modk, string base, Matriz2d* u){

    // * Função que gera um arquivo de dados no formato aceito pelo GNUPlot

    ofstream myfile;

    cout << "Gerando arquivo data" << to_string(k/modk) + base << "..." << endl;

    myfile.open("../data" + to_string(k/modk) + base);
        for (int j = STENCIL; j < d.Nz - STENCIL; j++){
            for (int i = STENCIL; i < d.Nx - STENCIL; i++){
                myfile << i*d.dz << " " << j*d.dx << " " << u->get(i, j) << "\n";
            }
            myfile << "\n\n";
        }
    myfile.close();
}

void salvaVTI(Dominio d, int k, int modk, Matriz2d* u){

    // * Função que gera um arquivo vtk ImageData para o ParaView

    ofstream myfile;

    cout << "Gerando arquivo data" << to_string(k/modk) << ".vti" << "..." << endl;

    myfile.open("../data" + to_string(k/modk) + ".vti");

    if(myfile.is_open()){

        myfile << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
        //                                          x1            x2          y1            y2          z1          z2
        myfile << "  <ImageData WholeExtent= \"" <<  0 << " " << d.Nx - 1 << " " << 0 << " " << d.Nz - 1 << " " << 0 << " " << 0 << "\" ";
        myfile << "Origin = \"" << 0 << " " << 0 << " " << 0 << "\" ";
        myfile << "Spacing = \"" << d.dx << " " << d.dz << " " << 0 << "\">\n";
        myfile << "    <Piece Extent = \"" << 0 << " " << d.Nx - 1 << " " << 0 << " " << d.Nz - 1 << " " << 0 << " " << 0 << "\">\n";
        myfile << "      <PointData Scalars=\"Amplitude\">\n";
        myfile << "        <DataArray type=\"Float32\" Name=\"Amplitude\" format=\"ascii\">\n";
        for (int j = 0; j < d.Nx; j++){
            for (int i = 0; i < d.Nz; i++){
                myfile << u->get(i, j) << " ";
            }
        }
        myfile << "\n        </DataArray>";
        myfile << "\n      </PointData>";
        myfile << "\n    </Piece>";
        myfile << "\n  </ImageData>";
        myfile << "\n</VTKFile>";

    } else {
        cout << "Erro na gravação do arquivo data" << to_string(k/modk) << ".vti" << endl;
    }

}

int main() {

    // inicializa os parametros da simulação
    Dominio d;
    leParametros(&d);

    // é necessária a utilização de apenas duas matrizes para implementar o método de diferenças finitas
    Matriz2d u_current(d.Nx, d.Nz);
    Matriz2d u_next(d.Nx, d.Nz);
    
    string base(".dat");
    int modk = 50;

    cout << "Nx = " << d.Nx << endl;
    cout << "Nz = " << d.Nz << endl;
    cout << "Nt = " << d.Nt << endl;
    cout << "Numero de Courant = " << d.cou << endl;

    auto inicio = chrono::high_resolution_clock::now();

    for (int k = 0; k <= d.Nt; k += 2){

        // calcula u_next
        mdfComAbsorcao(d, &u_current, &u_next, k);
        Reynolds(d, &u_current, &u_next);
        
        // calcula u_current
        mdfComAbsorcao(d, &u_next, &u_current, k + 1);
        // Reynolds(d, &u_next, &u_current);

        // * gera arquivo de dados a cada 100 iteracoes em k
        if (k % modk == 0){
            // salvaGNUPlot(d, k, modk, base, &u_current);
            salvaVTI(d, k, modk, &u_current);
        }

    }

    cout << "\nArquivos gerados com sucesso!" << endl;

    auto final = chrono::high_resolution_clock::now();
    chrono::duration<double> intervalo = final - inicio;
    cout << "\nTempo decorrido: " << intervalo.count() << "s\n";

    return 0;
}