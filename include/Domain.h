#include "Grid2d.h"
#include "Grid3d.h"

struct Domain {

    float X;       // * Distância em x em metros
    float Y;       // * Distância em y em metros
    float Z;       // * Distância em z em metros
    float T;       // * Tempo total da duração da simulação em segundos

    float c;       // * Celeridade da onda em m/s

    float dx;      // * Passo em x
    float dy;      // * Passo em y
    float dz;      // * Passo em z
    float dt;      // * Passo no tempo

    int Nx;        // * Numero de iteracoes em x
    int Ny;        // * Numero de iteracoes em y
    int Nz;        // * Numero de iteracoes em z
    int Nt;        // * Numero de iteracoes no tempo

    int n_fontes;  // * Quantidade de fontes
    float fcorte;  // * Frequencia de pico em Hz dos pulsos
    int xs;       // * posicao da fonte em x
    int ys;       // * posicao da fonte em y
    int zs;       // * posicao da fonte em z

    // * Velocidade variante
    Grid2d* vel;   // * Matriz 2D de velocidades

    // * Velocidade constante (obsoleto)
    float cou;     // * Número de courant, para dx = dz
    float c1; 
    float c2;

};
