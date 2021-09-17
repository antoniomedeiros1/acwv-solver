
#include <cstring>
using namespace std;

#include "../include/Grid2d.h"

Grid2d::Grid2d(int nl, int nc){
    this->n_linhas = nl;
    this->n_cols   = nc;
    int size       = this->n_linhas*this->n_cols;
    mat = new float[size];

    std::memset(mat,0.0,size*sizeof(float));

}

Grid2d::~Grid2d(){
    delete [] mat;
}

float Grid2d::get(int j, int i)
{
    return mat[j*this->n_cols + i];
}

void Grid2d::set(int j, int i, float val){
    mat[j*this->n_cols + i] = val;
}

float& Grid2d::operator()(int j, int i)
{
    return  mat[j*this->n_cols + i];
}

void Grid2d::swap(Grid2d &orig)
{
    int size = this->n_cols*this->n_linhas;
    for(int i = 0; i < size; i++)
    {
        float val    = this->mat[i];
        this->mat[i] = orig.mat[i];
        orig.mat[i]  = val;
    }
}
