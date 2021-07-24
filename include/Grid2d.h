// Use algumas extensao do vscode para ajudar na documentação de suas classes.
// Por exemplo: doxygen generator

// Pensou em usar classes templates?
// TODO: usar templates para utilizar diferentes tipos de dados 


class Grid2d{

    public:

        Grid2d(int n_linhas, int n_colunas);
        ~Grid2d();
        float get(int j, int i);
        void set(int j, int i, float val);
        void swap(Grid2d &orig);
        float& operator()(int j, int i); // utilizando sobrecarga de operadores

    private:
        int n_linhas;
        int n_cols;
        float* mat;

};