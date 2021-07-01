#define _USE_MATH_DEFINES
#define STENCIL 3 

#include <iostream>
#include <fstream>
#include <cmath>
#include "math.h"
#include <string>

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

float fonte(int x, int z, float t, float fcorte, float xs, float zs){

    // *funcao que simula um pulso sismico na posicao (xs, zs)

    float td = t - ((2*sqrt(M_PI))/fcorte);  
    float fc = (fcorte/(3*sqrt(M_PI)));

    if (x != xs || z != zs){
        return 0;
    } 

    float eq1 = (1.0 - 2.0 * M_PI * pow(M_PI * fc * td, 2));
    float eq2 = pow(M_E, M_PI*pow((M_PI*fc*td), 2));

    return eq1/eq2;

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
            + 2*u_current->get(i, j) - u_next->get(i, j) - d.c2 * fonte(i, j, k*d.dt, d.fcorte, d.xs/d.dx, d.zs/d.dz);

            u_next->set(i,j, val);
            
        }
    }

}

float atenuacao(float d){

    // * Função aplicada nas bordas para reduzir a amplitude da onda

    return pow(M_E, -1*pow(0.15*d, 2));

}

int dist(int i, int j, int Nx, int Nz){

    // * Essa funcao retorna a menor distancia em pontos entre o ponto (i, j) e o contorno 

    int minX, minZ;

    minX = min(i, Nx - i);
    minZ = min(j, Nz - j);

    return min(minX, minZ);

}

void geraArqGnuplot(Matriz2d u, int Nx, int Nz, int k, int modk, string base){

    // *essa funcao gera um arquivo de dados para plotar pelo GnuPlot
    // ! erro: core dump

    ofstream myfile;
    myfile.open("data" + to_string(k/modk) + base);
        for (int j = 0; j < Nz; j++){
            for (int i = 0; i < Nx; i++){
                myfile << i << " " << j << " " << u.get(i, j) << "\n";
            }
            myfile << "\n\n";
        }
    myfile.close();
}


int main() {

    Dominio d;
    leParametros(&d);

    float val; // variavel auxiliar 

    Matriz2d u_current(d.Nx, d.Nz);
    Matriz2d u_next(d.Nx, d.Nz);

    ofstream myfile;    
    string base(".dat");

    cout << "Nx = " << d.Nx << endl;
    cout << "Nz = " << d.Nz << endl;
    cout << "Nt = " << d.Nt << endl;
    cout << "Numero de Courant = " << d.cou << endl;

    for (int k = 0; k <= d.Nt; k += 2){

        // calcula u_next
        mdf(d, &u_current, &u_next, k);

        // * Reynolds NR BC

        // * borda esquerda
        //   du/dt - vel*du/dx = 0
        // (u(x,t+dt) - u(x,t))/dt = vel*(u(x+dx,t) - u(x,t))/dx
        // (u(x,t+dt) - u(x,t)) =   dt*(vel*(u(x+dx,t) - u(x,t))/dx
        // u(x, t+dt) = u(x,t) + cou * (u(x+dx,t) - u(x,t) )
        //    onde cou = dt*vel/dx
        for(int j = 0; j < d.Nz; j++) {
            u_next(STENCIL,j) = u_current(STENCIL,j) + d.cou*(u_current(STENCIL+1,j) - u_current(STENCIL,j));
        }

        // * borda direita
        //   du/dt + vel*du/dx = 0
        // u(x, t+dt) = u(x,t) - cou * (u(x,t) - u(x-dt,t) )
        // Aqui usamos diferencas atrasadas para discretizar do espaco
        for(int j = 0; j < d.Nz; j++) {
            u_next(d.Nx-STENCIL,j) = u_current(d.Nx-STENCIL,j) - d.cou*(u_current(d.Nx-STENCIL,j) - u_current(d.Nx-STENCIL-1,j));
        }

        // * borda superior
        //   du/dt - vel*du/dz = 0
        // u(z, t+dt) = u(z,t) + cou * (u(z+dz,t) - u(z,t) )
        for (int i = 0; i < d.Nx; i++) {
            u_next(i,STENCIL) = u_current(i,STENCIL) + d.cou*(u_current(i,STENCIL+1) - u_current(i,STENCIL));
            
        }

        // * borda inferior
        //   du/dt + vel*du/dz = 0
        // u(z, t+dt) = u(z,t) - cou * (u(z,t) - u(z-dz,t) )
        for (int i = 0; i < d.Nx; i++) {
            u_next(i,d.Nz-STENCIL) = u_current(i,d.Nz-STENCIL) - d.cou*(u_current(i,d.Nz-STENCIL) - u_current(i,d.Nz-STENCIL - 1));
        }
        
        // calcula u_current
        mdf(d, &u_next, &u_current, k);

        // * gera arquivo de dados a cada 100 iteracoes em k
        if (k % 100 == 0){
            cout << "Gerando arquivo data" << to_string(k/100) + base << "..." << endl;
            myfile.open(".././data" + to_string(k/100) + base);
            for (int i = 0; i < d.Nx; i++){
                for (int j = 0; j < d.Nz; j++)
                {
                    myfile << i*d.dx << " " << j*d.dz << " " << u_current(i, j) << "\n";
                }
                myfile << "\n\n";
            }
            myfile.close();
        }

    }

    cout << "Arquivos gerados com sucesso!" << endl;

    return 0;
}