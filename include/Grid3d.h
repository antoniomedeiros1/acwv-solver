
#include <cstring>

using namespace std;

class Grid3d{

    public:

        Grid3d(){}
        Grid3d(int nl, int nc, int np);
        ~Grid3d();
        float get(int k, int j, int i);
        void set(int k, int j, int i, float val);
        int getSize();
        float* firstptr();
        float& operator()(int k, int j, int i); // utilizando sobrecarga de operadores

    private:
        int n_linhas;
        int n_cols;
        int n_profundidade;
        float* mat;
        int size;
        
};