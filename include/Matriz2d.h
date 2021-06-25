// Use algumas extensao do vscode para ajudar na documentação de suas classes.
// Por exemplo: doxygen generator

// Pensou em usar classes templates?
// TODO: usar templates para utilizar diferentes tipos de dados 


class Matriz2d{

    public:

        Matriz2d(int n_linhas, int n_colunas);
        ~Matriz2d();
        float get(int i, int j);
        void set(int i, int j, float val);
        void swap(Matriz2d &orig);

        // Veja como podemos usar sobrecarga de operadores.
        float& operator()(int i, int j);

    private:
        int n_linhas;
        int n_cols;
        float* mat;

};