#include <iostream>
#include <fstream>
#include <cmath>
#include "math.h"
#include <string>
#include <chrono>
// #include "base64.h"
#include "Solver2d.h"

#define STENCIL 3

class Solver3d {

    public:

        // métodos
        Solver3d();  // ctor
        ~Solver3d(); // dtor
        void salvaVTI(Dominio d, Grid3d* u, string nomeDoArq, string info);
        void salvaVTIBinary(Dominio d, Grid3d* u, string nomeDoArq, string info);
        void solve();
    
    private:

        // métodos
        void leParametros(string nome);
        float pulso3d(int x, int y, int z, int t);
        void mdf3d(Grid3d* u1, Grid3d* u2, int t); // 
        void mdf3d(Grid3d* u1, Grid3d* u2);

        // atributos
        Dominio d;          // * Parametros do modelo
        Grid3d* u_current;  // * Matriz do estado atual 
        Grid3d* u_next;     // * Matriz do próximo estado

};