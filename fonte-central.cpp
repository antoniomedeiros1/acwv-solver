#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <string>
#include "Matriz3D.cpp"

#define M_PI 3.14159265358979323846 
#define M_e 2.71828182845904523536

using namespace std;

int main() {

    float X = 3000; //Largura do dominio em m
    float Z = 3000; //Altura do dominio em m
    float T = 1;    //Tempo total em s
    float c = 2200;  //Celeridade da onda em m/s
    float fcorte = 40;  //frequencia de pico = 40Hz

    float alpha = 5;
    float beta = 5;

    float dx = 10;  //Passos em x
    float dz = dx;  //Passos em z
    float dt = 0.00025;   //Passos no tempo

    int Nx = X/dx; //Iteracoes em x
    int Nz = Z/dz; //Iteracoes em z
    int Nt = 3500;//T/dt; //Iteracoes no tempo

    int xs = 150; //posicao da fonte em x
    int zs = 150; //posicao da fonte em z
    float cou = c*dt/dx; //numero de courant, para dx = dz

    cout << "Nx = " << Nx << endl;
    cout << "Nz = " << Nz << endl;
    cout << "Nt = " << Nt << endl;
    cout << "Numero de Courant = " << cou << endl;    

    float val;

    float fonte(int x, int z, float t, float fcorte, float xs, float zs);

    bool yet = false;

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

    }

    //gerar arquivos de dados para plotagem
    for(int k = 0; k < Nt; k += 50){
        myfile.open("/home/antonio/IC/modelagem_acustica_bidimensional/data/data" + to_string(k/50) + base);
        for (int j = 0; j < Nz; j++){
            for (int i = 0; i < Nx; i++){
                myfile << i << " " << j << " " << std::setprecision(17) << u.get(i, j, k) << "\n";
                /*
                if (i == 150 && j == 13){
                    cout << "u(" << i << ", " << j << ", " << k << ") = " << u.get(i, j, k) << endl;
                }*/
            }
            myfile << "\n\n";
        }
        myfile.close();
    }

}

float fonte(int x, int z, float t, float fcorte, float xs, float zs){

    float td = t - ((2*sqrt(M_PI))/fcorte); //para prevenir reflexao no tempo, isto eh, em z
    float fc = (fcorte/(3*sqrt(M_PI)));

    if (x != xs || z != zs){
        return 0;
    } 

    float eq1 = (1.0 - 2.0 * M_PI * pow(M_PI * fc * td, 2));

    float eq2 = pow(M_e, M_PI*pow((M_PI*fc*td), 2));

    return eq1/eq2;

}
