#include "Grid3d.h"

Grid3d::Grid3d(int nl, int nc, int np){

    this->n_linhas = nl;
    this->n_cols = nc;
    this->n_profundidade = np;

    this->size = this->n_linhas * this->n_cols * this->n_profundidade;

    this->mat = new float[size];

    memset(mat, 0.0, this->size*sizeof(float));

}

Grid3d::~Grid3d() {
    delete [] mat;
}

float Grid3d::get(int k, int j, int i) {
    return mat[k * this->n_cols * this->n_linhas + j * this->n_cols + i];
}

void Grid3d::set(int k, int j, int i, float val) {
    mat[k * this->n_cols * this->n_linhas + j * this->n_cols + i] = val;
}

int Grid3d::getSize(){
    return this->size;
}

float* Grid3d::firstptr(){
    return this->mat;
}

float& Grid3d::operator()(int k, int j, int i) {
    return  mat[k * this->n_cols * this->n_linhas + j * this->n_cols + i];
}