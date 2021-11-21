#include <iostream>
#include <fstream>
#include <cmath>
#include "math.h"
#include <string>
#include <chrono>

#include "Dominio.h"

#define STENCIL 3

using namespace std;

class Solver2d{

    public:

        // metodos
        Solver2d(string nomeDoArquivo);
        ~Solver2d();
        void imprimeParametros();
        void salvaVTI(Dominio d, Grid2d* u, string nomeDoArq, string info);
        void salvaVTIbin(Dominio d, Grid2d* u, string nomeDoArq, string info);
        void salvaSismograma(Grid2d sis, Dominio d);
        void solve();
        void salvaEstado();
        void carregaEstado();

        // atributos
    
    private: 

        // metodos
        void leParametros(string nome);
        void leModelo(string nome);

        void mdf(Dominio d, Grid2d* u_current, Grid2d* u_next, int k);
        void mdf(Dominio d, Grid2d* u_current, Grid2d* u_next);
        float fonte(int x, int z, float k);
        void aplicaReynolds(Grid2d* u_current, Grid2d* u_next);
        float atenuacao(float x, int borda);
        void aplicaAmortecimento();

        // atributos
        Dominio d;          // estrutura com os parametros 2d
        Grid2d* u_current;  // matriz de dados do estado atual
        Grid2d* u_next;     // matriz de dados do próximo estado
        Grid2d* f; 

        Grid2d* sis;        // matriz para registrar os traços sísmicos
        int posReceptor;    // posicao dos receptores em z

        string nome;        //  
        int modk;           // intervalo entre cada frame registrado

    
};
