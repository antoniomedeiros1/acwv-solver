class Matriz2d {

    public:

        Matriz2d(int x, int z);
        ~Matriz2d();
        float get(int i, int j);
        void set(int i, int j, float val);

    private:

        int l;
        int a;
        float* mat;

};