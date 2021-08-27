
class Grid2d{

    public:

        Grid2d(){}
        Grid2d(int n_linhas, int n_colunas);
        ~Grid2d();
        float get(int j, int i);
        void set(int j, int i, float val);
        void swap(Grid2d &orig);
        float& operator()(int j, int i); // utilizando sobrecarga de operadores
        float* firstptr(){ return mat; }
        int getSize(){ return n_linhas * n_cols; };

    private:
        int n_linhas;
        int n_cols;
        float* mat;

};