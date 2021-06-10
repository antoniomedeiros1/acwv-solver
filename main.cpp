#define _USE_MATH_DEFINES
#define STENCIL 3 

#include <iostream>
#include <fstream>
#include <cmath>
#include "math.h"
#include <iomanip>
#include <string>
#include "Matriz3D.cpp"
#include "./include/Matriz2d.h"

using namespace std;

// 1. Melhorar a estrutura do Codigo
//   1.a: Utilize estrutura de diretorio do projeto:
//        src/            --> todos os arquivos .cpp
//        include/        --> todos os arquivos .h
//        main.cpp
//        CMakeLists.txt  --> Arquivo de configuração de projeto qye usa Cmake.
// 2. Estudar https://cmake.org/ e gerar um scrip CMakeLists.txt para compilar 
//       o código considerando a estrutura acima.
// 3. Modificar o código para trabalhar apenas com dois/tres vetores.
// 4. Considerar o tamanho do stencil nos loops i,j
// 5. Implementar a camada de atenuação por borda.
// 6: Remover os ifs internos nas funções get de matriz3D e atenuação
// 7. Criar arquivo com funções de escrita e chamar a função dentro do laco do tempo.
// 8. Implementar uma estrutura ou classe para armazenar os parametros do problema

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

    float fonte(int x, int z, float t, float fcorte, float xs, float zs);
    float atenuacao(int x, int z, int Nx, int Nz);//funcao usada para impedir reflexoes nas bordas
    void plot1d(int t1, int t2, int passo, Matriz3D u, int Nx);

    bool yet = false;

    Matriz2d u_current(Nx, Nz);
    Matriz2d u_old(Nx, Nz);
    Matriz2d u_next(Nx, Nz);

    for (int k = 0; k < Nt; k++){
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
                + 2*u_current.get(i, j) - u_next.get(i, j) - 
                c2 * fonte(i, j, k*dt, fcorte, xs, zs);

                u_next.set(i, j, val);

            }

        }

        // gera arquivo com os dados de u_next

        // atualiza
    }

    Matriz3D u(Nx, Nz, Nt);
    ofstream myfile;    
    string base(".dat");
    
    //MDF
    for (int k = -1; k < Nt - 1; k++){

        for (int j = 0; j < Nz; j++){

            for (int i = 0; i < Nx; i++){

                val = 
                (pow(cou, 2)/12.0) *
                (
                    -1*(u.get(i - 2, j, k) + u.get(i, j - 2, k)) + 
                    16*(u.get(i - 1, j, k) + u.get(i, j - 1, k)) - 
                    60 * u.get(i, j, k) +
                    16*(u.get(i + 1, j, k) + u.get(i, j + 1, k)) -
                    (u.get(i + 2, j, k) + u.get(i, j + 2, k)) 
                )
                + 2*u.get(i, j, k) - u.get(i, j, k - 1) - 
                pow(c*dt, 2) * fonte(i, j, k*dt, fcorte, xs, zs);

                u.set(i, j, k + 1, val);

            }

        }

        for (int j = 0; j < Nz; j++){
            for (int i = 0; i < Nx; i++) {
                u.set(i, j, k + 1, u.get(i, j, k + 1)*atenuacao(i, j, Nx, Nz));
            }
        }

    }

    //gerar arquivos de dados para plotagem
    for(int k = 0; k < Nt; k += 50){
        myfile.open("/home/antonio/IC/modelagem_acustica_bidimensional/data/data" + to_string(k/50) + base);
        for (int j = 0; j < Nz; j++){
            for (int i = 0; i < Nx; i++){
                myfile << i << " " << j << " " << setprecision(17) << u.get(i, j, k) << "\n";
                /*
                if (i == 150 && j == 13){
                    cout << "u(" << i << ", " << j << ", " << k << ") = " << u.get(i, j, k) << endl;
                }*/
            }
            myfile << "\n\n";
        }
        myfile.close();
    }

    //plot1d(2400,Nt, 5, u, Nx);
}

void plot1d(int t1, int t2, int passo, Matriz3D u, int Nx){

    ofstream myfile;

    myfile.open("/home/antonio/IC/modelagem_acustica_bidimensional/data1d/data.dat");
    for(int k = t1; k < t2; k += passo){

        for (int i = 0; i < Nx; i++){
            myfile << i << " " << std::setprecision(17) << u.get(i, 150, k) << "\n";
            
        }
        myfile << "\n\n";
    }
    myfile.close();

}

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
