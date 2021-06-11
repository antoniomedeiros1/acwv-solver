#define _USE_MATH_DEFINES
#define STENCIL 3 

#include <iostream>
#include <fstream>
#include <cmath>
#include "math.h"
#include <iomanip>
#include <string>
#include "./include/Matriz2d.h"

using namespace std;

// 1. Melhorar a estrutura do Codigo
//   1.a: Utilize estrutura de diretorio do projeto:
//        src/            --> todos os arquivos .cpp ✔
//        include/        --> todos os arquivos .h ✔
//        main.cpp ✔
//        CMakeLists.txt  --> Arquivo de configuração de projeto qye usa Cmake. ✔
// 2. Estudar https://cmake.org/ e gerar um scrip CMakeLists.txt para compilar ✔
//       o código considerando a estrutura acima.
// 3. Modificar o código para trabalhar apenas com dois/tres vetores. ✔
// 4. Considerar o tamanho do stencil nos loops i,j ✔
// 5. Implementar a camada de atenuação por borda.
// 6: Remover os ifs internos nas funções get de matriz3D e atenuação ✔
// 7. Criar arquivo com funções de escrita e chamar a função dentro do laco do tempo.
// 8. Implementar uma estrutura ou classe para armazenar os parametros do problema


float fonte(int x, int z, float t, float fcorte, float xs, float zs){

    float td = t - ((2*sqrt(M_PI))/fcorte);
    float fc = (fcorte/(3*sqrt(M_PI)));

    if (x != xs || z != zs){
        return 0;
    } 

    float eq1 = (1.0 - 2.0 * M_PI * pow(M_PI * fc * td, 2));

    float eq2 = pow(M_E, M_PI*pow((M_PI*fc*td), 2));

    return eq1/eq2;

}

float atenuacao(int x, int z, int Nx, int Nz){

    int n = 20;
    if (x < n){
        return pow(M_E, -1*pow(0.098*(n - x), 2));
    } else if (x > Nx - n){
        return pow(M_E, -1*pow(0.098*(n - (Nx - x)), 2));
    } else if (z < n){
        return pow(M_E, -1*pow(0.098*(n - z), 2));
    } else if (z > Nz - n){
        return pow(M_E, -1*pow(0.098*(n - (Nz - z)), 2));
    } else {
        return 1.0;
    }

}

int main() {

    float X = 3000; //Largura do dominio em m
    float Z = 3000; //Altura do dominio em m
    float T = 1;    //Tempo total em s
    float c = 2200;  //Celeridade da onda em m/s
    float fcorte = 40;  //frequencia de pico = 40Hz
    float dx = 10;  //Passos em x
    float dz = dx;  //Passos em z
    float dt = 0.00025;   //Passos no tempo

    int Nx = X/dx; //Iteracoes em x
    int Nz = Z/dz; //Iteracoes em z
    int Nt = 3500;//T/dt; //Iteracoes no tempo

    int xs = 150; //posicao da fonte em x
    int zs = 150; //posicao da fonte em z
    float cou = c*dt/dx; //numero de courant, para dx = dz
    float c1 = (pow(cou, 2)/12.0);
    float c2 = pow(c*dt, 2);

    cout << "Nx = " << Nx << endl;
    cout << "Nz = " << Nz << endl;
    cout << "Nt = " << Nt << endl;
    cout << "Numero de Courant = " << cou << endl;    

    float val;
    float w; //fator atenuador

    Matriz2d u_current(Nx, Nz);
    Matriz2d u_old(Nx, Nz);
    Matriz2d u_next(Nx, Nz);

    ofstream myfile;    
    string base(".dat");

    for (int k = 0; k < /*Nt*/ 10; k++){

        // calcula u_next
        for (int j = STENCIL; j < Nz - STENCIL; j++){
            for (int i = STENCIL; i < Nx - STENCIL; i++){

                val = 
                c1 *
                (
                    -1*(u_current.get(i - 2, j) + u_current.get(i, j - 2)) + 
                    16*(u_current.get(i - 1, j) + u_current.get(i, j - 1)) - 
                    60* u_current.get(i, j) +
                    16*(u_current.get(i + 1, j) + u_current.get(i, j + 1)) -
                    (u_current.get(i + 2, j) + u_current.get(i, j + 2)) 
                ) 
                + 2*u_current.get(i, j) - u_next.get(i, j) - c2 * fonte(i, j, k*dt, fcorte, xs, zs);

                u_next.set(i, j, val);
            }
        }

        // troca u_next <--> u_current
        for (int j = STENCIL; j < Nz - STENCIL; j++){
            for (int i = STENCIL; i < Nx - STENCIL; i++){
                val = u_next.get(i, j);
                u_next.set(i, j, u_current.get(i, j));
                u_current.set(i, j, val);
            }
        }

        // gera arquivo de dados a cada 50 iteracoes em k
        if (k % 5 == 0){
            myfile.open(".././data" + to_string(k/5) + base);
            for (int j = 0; j < Nz; j++){
                for (int i = 0; i < Nx; i++){
                    myfile << i << " " << j << " " << setprecision(17) << u_current.get(i, j) << "\n";
                }
                myfile << "\n\n";
            }
            myfile.close();
        }

    }
}

