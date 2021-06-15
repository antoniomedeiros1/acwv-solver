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
// 5. Implementar a camada de atenuação por borda. (implementado porem nao esta funcionando corretamente)
// 6: Remover os ifs internos nas funções get de matriz3D e atenuação ✔
// 7. Criar arquivo com funções de escrita e chamar a função dentro do laco do tempo. ✔
// 8. Implementar uma estrutura ou classe para armazenar os parametros do problema


float fonte(int x, int z, float t, float fcorte, float xs, float zs){

    // *funcao que simula um pulso sismico na posica (xs, zs)

    float td = t - ((2*sqrt(M_PI))/fcorte);  
    float fc = (fcorte/(3*sqrt(M_PI)));

    if (x != xs || z != zs){
        return 0;
    } 

    float eq1 = (1.0 - 2.0 * M_PI * pow(M_PI * fc * td, 2));

    float eq2 = pow(M_E, M_PI*pow((M_PI*fc*td), 2));

    return eq1/eq2;

}

float atenuacao(int d){

    // *funcao que recebe a distancia entre um ponto na borda e o contorno e retorna um termo atenuador

    return pow(M_E, -1*pow(0.098*d, 2));
}

void geraArqGnuplot(Matriz2d u, int Nx, int Nz, int k, int modk, string base){

    // *essa funcao gera um arquivo de dados para plotar pelo GnuPlot
    // ! erro: core dump

    ofstream myfile;
    myfile.open("./data" + to_string(k/modk) + base);
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
    int Nt = 3500;      // T/dt; //Iteracoes no tempo

    int xs = 150;       // posicao da fonte em x
    int zs = 150;       // posicao da fonte em z

    float cou = c*dt/dx;  // numero de courant, para dx = dz
    float c1 = (pow(cou, 2)/12.0); 
    float c2 = pow(c*dt, 2);
    int borda = 20;       // largura da borda que sofrera a atenuacao
    int d;                // distancia entre o no e o contorno
    float atenuacoes[borda];  

    float val; // variavel auxiliar 

    // para calcular o mdf basta utilizar duas matrizes, uma para armazenar o valor atual
    // e outra para calcular o proximo valor 
    // como essa forma discreta da equacao da onda requer o valor dois passos atras no tempo (i, j, k - 2)
    // basta realizar a troca dos valores de u_current para a u_next, assim na proxima iteracao os valores
    // no tempo k - 2 ficam armazenados na mesma matriz que armazenara os valores futuros
    Matriz2d u_current(Nx, Nz);
    Matriz2d u_next(Nx, Nz);

    ofstream myfile;    
    string base(".dat");

    // *imprime os parametros 
    cout << "Nx = " << Nx << endl;
    cout << "Nz = " << Nz << endl;
    cout << "Nt = " << Nt << endl;
    cout << "Numero de Courant = " << cou << endl;

    // *calcula o vetor de atenuacoes para evitar que o calculo seja realizada a cada iteracao na borda
    for(int i = 0; i < borda; i++){
        atenuacoes[i] = atenuacao(i);
    }

    // *MDF 
    // *o algoritmo segue a seguinte ordem:
    // *     1. inicializa duas matrizes zeradas u_current e u_next
    // *     2. calcula u_next a partir dos valores da u_current pelo MDF
    // *     3. percorre as bordas aplicando uma funcao de atenuacao do valor para evitar a reflexao
    // *     4. u_current passa a ser u_next e u_next recebe os valores anteriores de u_current
    // *     5. gera um arquivo de dados para plot a cada x iteracoes no tempo
    for (int k = 0; k < Nt; k++){

        // *calcula u_next
        for (int j = STENCIL; j < Nz - STENCIL; j++){
            for (int i = STENCIL; i < Nx - STENCIL; i++){

                val = 
                c1 *
                (
                    -1*(u_current.get(i - 2, j) + u_current.get(i, j - 2)) + 
                    16*(u_current.get(i - 1, j) + u_current.get(i, j - 1)) - 
                    60* u_current.get(i,     j) +
                    16*(u_current.get(i + 1, j) + u_current.get(i, j + 1)) -
                       (u_current.get(i + 2, j) + u_current.get(i, j + 2)) 
                ) 
                + 2*u_current.get(i, j) - u_next.get(i, j) - c2 * fonte(i, j, k*dt, fcorte, xs, zs);

                u_next.set(i, j, val);
            }
        }

        // *percorre as bordas do modelo aplicando a atenuacao
        for (int j = 0; j < borda; j++){
            for(int i = 0; i < Nx; i++){
                u_next.set(i, j, u_next.get(i, j) * atenuacoes[borda - j]);
            }
        }

        // *troca u_next <--> u_current
        for (int j = STENCIL; j < Nz - STENCIL; j++){
            for (int i = STENCIL; i < Nx - STENCIL; i++){
                val = u_next.get(i, j);
                u_next.set(i, j, u_current.get(i, j));
                u_current.set(i, j, val);
            }
        }

        // *gera arquivo de dados a cada 500 iteracoes em k
        if (k % 200 == 0){
            myfile.open(".././data" + to_string(k/200) + base);
            for (int j = 0; j < Nz; j++){
                for (int i = 0; i < Nx; i++){
                    myfile << i << " " << j << " " << u_current.get(i, j) << "\n";
                }
                myfile << "\n\n";
            }
            myfile.close();
        }

    }
}

