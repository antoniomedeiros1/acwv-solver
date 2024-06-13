#include <petscdm.h>
#include <petscdmda.h>

// struct Domain {
//     float X;
//     float Y;
//     float Z;
//     float T;
//     float c;
//     float dx;
//     float dy;
//     float dz;
//     float dt;
//     int Nx;
//     int Ny;
//     int Nz;
//     int Nt;
//     int n_fontes;
//     float fcorte;
//     int xs;
//     int ys;
//     int zs;
//     Grid2d* vel;
//     float cou;
//     float c1; 
//     float c2;
// };

struct Domain {
    PetscReal X, Y, Z, T, c, dx, dy, dz, dt;
    PetscInt Nx, Ny, Nz, Nt, xs, ys, zs;
    Vec vel;
    Vec vel_local;
    PetscReal cou, c1, c2;
};
