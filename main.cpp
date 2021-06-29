#define _USE_MATH_DEFINES
#define STENCIL 3 

#include <iostream>
#include <fstream>
#include <cmath>
#include "math.h"
#include <iomanip>
#include <string>

#include "Matriz2d.h"

using namespace std;

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

    float X = 3000;     // Largura do dominio em m
    float Z = 3000;     // Altura do dominio em m
    float T = 1;        // Tempo total em s
    float c = 2200;     // Celeridade da onda em m/s
    float fcorte = 40;  // frequencia de pico em Hz
    float dx = 10;      // Passos em x
    float dz = dx;      // Passos em z
    float dt = 0.00025; // Passos no tempo

    int Nx = X/dx;      // Iteracoes em x
    int Nz = Z/dz;      // Iteracoes em z
    int Nt = T/dt;      //Iteracoes no tempo

    int xs = 150;       // posicao da fonte em x
    int zs = 150;       // posicao da fonte em z

    float cou = c*dt/dx;  // numero de courant, para dx = dz
    float c1 = (pow(cou, 2)/12.0); 
    float c2 = pow(c*dt, 2);

    // int borda = 20;       // largura da borda que sofrera a atenuacao

    float val; // variavel auxiliar 

    Matriz2d u_current(Nx, Nz);
    Matriz2d u_next(Nx, Nz);

    ofstream myfile;    
    string base(".dat");

    // *imprime os parametros 
    cout << "Nx = " << Nx << endl;
    cout << "Nz = " << Nz << endl;
    cout << "Nt = " << Nt << endl;
    cout << "Numero de Courant = " << cou << endl;

    // * MDF 
    // * o algoritmo segue a seguinte ordem:
    // *     1. inicializa duas matrizes zeradas u_current e u_next
    // *     2. calcula u_next a partir dos valores da u_current pelo MDF
    // *     3. percorre as bordas aplicando uma funcao de atenuacao do valor para evitar a reflexao
    // *     4. u_current passa a ser u_next e u_next recebe os valores anteriores de u_current
    // *     5. gera um arquivo de dados para plot a cada x iteracoes no tempo
    for (int k = 0; k <= Nt; k++){

        // calcula u_next
        for (int i = STENCIL; i <= Nx - STENCIL; i++){
            for (int j = STENCIL; j <= Nz - STENCIL; j++){

                // if (dist(i, j, Nx, Nz) > borda){

                // *calcula o MDF normalmente

                val = 
                c1 *
                (
                    -1*(u_current(i - 2, j) + u_current(i, j - 2)) + 
                    16*(u_current(i - 1, j) + u_current(i, j - 1)) - 
                    60* u_current(i,     j) +
                    16*(u_current(i + 1, j) + u_current(i, j + 1)) -
                    (u_current(i + 2, j) + u_current(i, j + 2)) 
                ) 
                + 2*u_current(i, j) - u_next(i, j) - c2 * fonte(i, j, k*dt, fcorte, xs, zs);

                u_next(i,j) = val;

                // } else {

                //     // *calcula o MDF e aplica uma funcao de atenuacao 
                //     val = 
                //     (c1 *
                //     (
                //         -1*(u_current(i - 2, j) + u_current(i, j - 2)) + 
                //         16*(u_current(i - 1, j) + u_current(i, j - 1)) - 
                //         60* u_current(i,     j) +
                //         16*(u_current(i + 1, j) + u_current(i, j + 1)) -
                //         (u_current(i + 2, j) + u_current(i, j + 2)) 
                //     ) 
                //     + 2*u_current(i, j) - u_next(i, j));

                //     val *= atenuacao(dist(i, j, Nx, Nz));

                //     u_next(i,j) = val;

                // }
                
            }
        }

        // * Reynolds NR BC

        // * borda esquerda
        //   du/dt - vel*du/dx = 0
        // (u(x,t+dt) - u(x,t))/dt = vel*(u(x+dx,t) - u(x,t))/dx
        // (u(x,t+dt) - u(x,t)) =   dt*(vel*(u(x+dx,t) - u(x,t))/dx
        // u(x, t+dt) = u(x,t) + cou * (u(x+dx,t) - u(x,t) )
        //    onde cou = dt*vel/dx
        for(int j = 0; j < Nz; j++) {
            u_next(STENCIL,j) = u_current(STENCIL,j) + cou*(u_current(STENCIL+1,j) - u_current(STENCIL,j));
        }

        // * borda direita
        //   du/dt + vel*du/dx = 0
        // u(x, t+dt) = u(x,t) - cou * (u(x,t) - u(x-dt,t) )
        // Aqui usamos diferencas atrasadas para discretizar do espaco
        for(int j = 0; j < Nz; j++) {
            u_next(Nx-STENCIL,j) = u_current(Nx-STENCIL,j) - cou*(u_current(Nx-STENCIL,j) - u_current(Nx-STENCIL-1,j));
        }

        // * borda superior
        //   du/dt - vel*du/dz = 0
        // u(z, t+dt) = u(z,t) + cou * (u(z+dz,t) - u(z,t) )
        for (int i = 0; i < Nx; i++) {
            u_next(i,STENCIL) = u_current(i,STENCIL) + cou*(u_current(i,STENCIL+1) - u_current(i,STENCIL));
            
        }

        // * borda inferior
        //   du/dt + vel*du/dz = 0
        // u(z, t+dt) = u(z,t) - cou * (u(z,t) - u(z-dz,t) )
        for (int i = 0; i < Nx; i++) {
            u_next(i,Nz-STENCIL) = u_current(i,Nz-STENCIL) - cou*(u_current(i,Nz-STENCIL) - u_current(i,Nz-STENCIL - 1));
        }
        
        // TODO: Essa troca de vetores pode ser evitada
        u_next.swap(u_current);

        // * gera arquivo de dados a cada 100 iteracoes em k
        if (k % 100 == 0){
            cout << "Gerando arquivo data" << to_string(k/100) + base << "..." << endl;
            myfile.open(".././data" + to_string(k/100) + base);
            for (int i = 0; i < Nx; i++){
                for (int j = 0; j < Nz; j++)
                {
                    myfile << i*dx << " " << j*dz << " " << u_current(i, j) << "\n";
                }
                myfile << "\n\n";
            }
            myfile.close();
        }

    }

    cout << "Arquivos gerados com sucesso!" << endl;

    return 0;
}