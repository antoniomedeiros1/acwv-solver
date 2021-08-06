struct Dominio {

    float X;       // *Largura do dominio em m
    float Z;       // *Altura do dominio em m
    float T;       // *Tempo total em s

    float c;       // *Celeridade da onda em m/s

    float dx;      // *Passos em x
    float dz;      // *Passos em z
    float dt;      // *Passos no tempo

    int Nx;        // *Numero de iteracoes em x
    int Nz;        // *Numero de iteracoes em z
    int Nt;        // *T/dt; Numero de iteracoes no tempo

    float fcorte;  // *Frequencia de pico em Hz
    int xs;     // *posicao da fonte em x
    int zs;     // *posicao da fonte em z

    Grid2d* vel;

    float cou;  // numero de courant, para dx = dz
    float c1; 
    float c2;

};
