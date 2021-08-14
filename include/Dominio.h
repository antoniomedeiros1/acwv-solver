#include "Grid2d.h"
#include "Grid3d.h"
struct Dominio {

    float X;       // *Distância em x em metros
    float Y;       // *Distância em y em metros
    float Z;       // *Distância em z em metros
    float T;       // *Tempo total da duração da simulação em segundos

    float c;       // *Celeridade da onda em m/s

    float dx;      // *Passo em x
    float dy;      // *Passo em y
    float dz;      // *Passo em z
    float dt;      // *Passo no tempo

    int Nx;        // *Numero de iteracoes em x
    int Ny;        // *Numero de iteracoes em y
    int Nz;        // *Numero de iteracoes em z
    int Nt;        // *T/dt; Numero de iteracoes no tempo

    float fcorte;  // *Frequencia de pico em Hz
    int xs;        // *posicao da fonte em x
    int ys;        // *posicao da fonte em y
    int zs;        // *posicao da fonte em z

    // *Velocidade variavel
    Grid2d* vel;   // *Matriz 2D de velocidades

    // *Velocidade constante
    float cou;     // *numero de courant, para dx = dz
    float c1; 
    float c2;

};
