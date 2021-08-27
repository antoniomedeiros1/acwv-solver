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
        void solve();

        // atributos
    
    private: 

        // metodos
        void leParametros(string nome);

        void mdf(Dominio d, Grid2d* u_current, Grid2d* u_next, int k);
        float fonte(int x, int z, float k);
        void aplicaReynolds();
        float atenuacao(float x, int borda);
        void aplicaAmortecimento();

        // atributos
        Dominio d; // estrutura com os parametros 2d
        Grid2d* u_current;
        Grid2d* u_next;

        Grid2d* sis; // matriz para registrar os traços sísmicos
        int posReceptor; // posicao dos receptores em z

        string nome;
        int modk; // intervalo entre cada frame registrado

    
};