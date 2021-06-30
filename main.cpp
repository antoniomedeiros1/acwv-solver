#define _USE_MATH_DEFINES
#define STENCIL 3 

#include <iostream>
#include <fstream>
#include <cmath>
#include "math.h"
#include <iomanip>
#include <string>

#include "Matriz2d.h"
#include "Domain.h"

using namespace std;

void setParameters(Domain* dom){

    ifstream myfile;
    
    myfile.open(".././parametros.txt", ios::in);

    if (myfile.is_open()){

        cout << "Lendo parametros..." << endl;

        myfile >> dom->X;       // *Largura do dominio em m
        myfile >> dom->Z;       // *Altura do dominio em m
        myfile >> dom->T;       // *Tempo total em s
        myfile >> dom->c;       // *Celeridade da onda em m/s
        myfile >> dom->fcorte;  // *Frequencia de pico em Hz
        myfile >> dom->dx;      // *Passos em x
        dom->dz = dom->dx;      // *Passos em z
        myfile >> dom->dt;      // *Passos no tempo

        dom->Nx = dom->X/dom->dx;        // *Numero de iteracoes em x
        dom->Nz = dom->Z/dom->dz;        // *Numero de iteracoes em z
        dom->Nt = dom->T/dom->dt;        // *T/dt; Numero de iteracoes no tempo

        myfile >> dom->xs;     // *posicao da fonte em x
        myfile >> dom->zs;     // *posicao da fonte em z
        
        dom->cou = dom->c*dom->dt/dom->dx;  // numero de courant, para dx = dz
        dom->c1 = (pow(dom->cou, 2)/12.0); 
        dom->c2 = pow(dom->c*dom->dt, 2);

    } else {
        cout << "Falha ao abrir arquivo de parametros" << endl;
    }

    myfile.close();

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

float atenuacao(float d){
    return pow(M_E, -1*pow(0.15*d, 2));
}

int dist(int i, int j, int Nx, int Nz){

    // * Essa funcao retorna a menor distancia em pontos entre o ponto (i, j) e as bordas 

    int minX, minZ;

    minX = min(i, Nx - i);
    minZ = min(j, Nz - j);

    return min(minX, minZ);

}

void geraArqGnuplot(Matriz2d u, int Nx, int Nz, int k, int modk, string base){

    // *essa funcao gera um arquivo de dados para plotar pelo GnuPlot
    // ! erro: core dump

    ofstream myfile;
    myfile.open(".././data" + to_string(k/modk) + base);
        for (int j = 0; j < Nz; j++){
            for (int i = 0; i < Nx; i++){
                myfile << i << " " << j << " " << u.get(i, j) << "\n";
            }
            myfile << "\n\n";
        }
    myfile.close();
}


int main() {

    Domain dom;

    float val; // variavel auxiliar 

    Matriz2d u_current(dom.Nx, dom.Nz);
    Matriz2d u_next(dom.Nx, dom.Nz);

    ofstream myfile;    
    string base(".dat");

    setParameters(&dom);

    // * imprime os parametros 
    cout << "Nx = " << dom.Nx << endl;
    cout << "Nz = " << dom.Nz << endl;
    cout << "Nt = " << dom.Nt << endl;
    cout << "Numero de Courant = " << dom.cou << endl;

    for (int k = 0; k <= dom.Nt; k++){

        // calcula u_next
        for (int i = STENCIL; i <= dom.Nx - STENCIL; i++){
            for (int j = STENCIL; j <= dom.Nz - STENCIL; j++){

                val = 
                dom.c1 *
                (
                    -1*(u_current(i - 2, j) + u_current(i, j - 2)) + 
                    16*(u_current(i - 1, j) + u_current(i, j - 1)) - 
                    60* u_current(i,     j) +
                    16*(u_current(i + 1, j) + u_current(i, j + 1)) -
                    (u_current(i + 2, j) + u_current(i, j + 2)) 
                ) 
                + 2*u_current(i, j) - u_next(i, j) - dom.c2 * fonte(i, j, k*dom.dt, dom.fcorte, dom.xs, dom.zs);

                u_next(i,j) = val;
                
            }
        }

        // * Reynolds NR BC

        // * borda esquerda
        //   du/dt - vel*du/dx = 0
        // (u(x,t+dt) - u(x,t))/dt = vel*(u(x+dx,t) - u(x,t))/dx
        // (u(x,t+dt) - u(x,t)) =   dt*(vel*(u(x+dx,t) - u(x,t))/dx
        // u(x, t+dt) = u(x,t) + cou * (u(x+dx,t) - u(x,t) )
        //    onde cou = dt*vel/dx
        for(int j = 0; j < dom.Nz; j++) {
            u_next(STENCIL,j) = u_current(STENCIL,j) + dom.cou*(u_current(STENCIL+1,j) - u_current(STENCIL,j));
        }

        // * borda direita
        //   du/dt + vel*du/dx = 0
        // u(x, t+dt) = u(x,t) - cou * (u(x,t) - u(x-dt,t) )
        // Aqui usamos diferencas atrasadas para discretizar do espaco
        for(int j = 0; j < dom.Nz; j++) {
            u_next(dom.Nx-STENCIL,j) = u_current(dom.Nx-STENCIL,j) - dom.cou*(u_current(dom.Nx-STENCIL,j) - u_current(dom.Nx-STENCIL-1,j));
        }

        // * borda superior
        //   du/dt - vel*du/dz = 0
        // u(z, t+dt) = u(z,t) + cou * (u(z+dz,t) - u(z,t) )
        for (int i = 0; i < dom.Nx; i++) {
            u_next(i,STENCIL) = u_current(i,STENCIL) + dom.cou*(u_current(i,STENCIL+1) - u_current(i,STENCIL));
            
        }

        // * borda inferior
        //   du/dt + vel*du/dz = 0
        // u(z, t+dt) = u(z,t) - cou * (u(z,t) - u(z-dz,t) )
        for (int i = 0; i < dom.Nx; i++) {
            u_next(i,dom.Nz-STENCIL) = u_current(i,dom.Nz-STENCIL) - dom.cou*(u_current(i,dom.Nz-STENCIL) - u_current(i,dom.Nz-STENCIL - 1));
        }
        
        // TODO: Essa troca de vetores pode ser evitada
        u_next.swap(u_current);

        // * gera arquivo de dados a cada 100 iteracoes em k
        if (k % 100 == 0){
            cout << "Gerando arquivo data" << to_string(k/100) + base << "..." << endl;
            myfile.open(".././data" + to_string(k/100) + base);
            for (int i = 0; i < dom.Nx; i++){
                for (int j = 0; j < dom.Nz; j++)
                {
                    myfile << i*dom.dx << " " << j*dom.dz << " " << u_current(i, j) << "\n";
                }
                myfile << "\n\n";
            }
            myfile.close();
        }

    }

    cout << "Arquivos gerados com sucesso!" << endl;

    return 0;
}