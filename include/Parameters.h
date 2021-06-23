
class Parameters {

    public:

        Parameters();
        ~Parameters();

    private:

        float X;       // *Largura do dominio em m
        float Z;       // *Altura do dominio em m
        float T;       // *Tempo total em s
        float c;       // *Celeridade da onda em m/s
        float fcorte;  // *Frequencia de pico em Hz
        float dx;      // *Passos em x
        float dz;      // *Passos em z
        float dt;      // *Passos no tempo
        int Nx;        // *Numreo de iteracoes em x
        int Nz;        // *Numero de iteracoes em z
        int Nt;        // *T/dt; Numero de iteracoes no tempo

        int xs;     // *posicao da fonte em x
        int zs;     // *posicao da fonte em z

        float cou;  // *numero de courant, para dx = dz
        float c1;   // *constante (1) utilizada no mdf 
        float c2;   // *constante (2) utilizada no mdf 
        int borda;  // *largura da borda que sofrera a atenuacao

}