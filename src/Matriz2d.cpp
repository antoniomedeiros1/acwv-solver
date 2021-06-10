#include "../include/Matriz2d.h"

Matriz2d::Matriz2d(int x, int z){
    l = x;
    a = z;
    mat = new float[l*a];
    for (int i = 0; i < l*a; i++){
        mat[i] = 0;
    }
}

Matriz2d::~Matriz2d(){
    delete mat;
}

float Matriz2d::get(int i, int j){
    return mat[j*l + i];
}

void Matriz2d::set(int i, int j, float val){
    mat[j*l + i] = val;
}
