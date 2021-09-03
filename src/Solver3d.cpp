#include "Solver3d.h"


const static char* b64="ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/" ;

// Converts binary data of length=len to base64 characters.
// Length of the resultant string is stored in flen
// (you must pass pointer flen).

char* base64( const void* binaryData, int len, int *flen )
{
  const unsigned char* bin = (const unsigned char*) binaryData ;
  char* res ;
  
  int rc = 0 ; // result counter
  int byteNo ; // I need this after the loop
  
  int modulusLen = len % 3 ;
  int pad = ((modulusLen&1)<<1) + ((modulusLen&2)>>1) ; // 2 gives 1 and 1 gives 2, but 0 gives 0.
  
  *flen = 4*(len + pad)/3 ;
  res = (char*) malloc( *flen + 1 ) ; // and one for the null
  if( !res )
  {
    puts( "ERROR: base64 could not allocate enough memory." ) ;
    puts( "I must stop because I could not get enough" ) ;
    return 0;
  }
  
  for( byteNo = 0 ; byteNo <= len-3 ; byteNo+=3 )
  {
    unsigned char BYTE0=bin[byteNo];
    unsigned char BYTE1=bin[byteNo+1];
    unsigned char BYTE2=bin[byteNo+2];
    res[rc++]  = b64[ BYTE0 >> 2 ] ;
    res[rc++]  = b64[ ((0x3&BYTE0)<<4) + (BYTE1 >> 4) ] ;
    res[rc++]  = b64[ ((0x0f&BYTE1)<<2) + (BYTE2>>6) ] ;
    res[rc++]  = b64[ 0x3f&BYTE2 ] ;
  }
  
  if( pad==2 )
  {
    res[rc++] = b64[ bin[byteNo] >> 2 ] ;
    res[rc++] = b64[ (0x3&bin[byteNo])<<4 ] ;
    res[rc++] = '=';
    res[rc++] = '=';
  }
  else if( pad==1 )
  {
    res[rc++]  = b64[ bin[byteNo] >> 2 ] ;
    res[rc++]  = b64[ ((0x3&bin[byteNo])<<4)   +   (bin[byteNo+1] >> 4) ] ;
    res[rc++]  = b64[ (0x0f&bin[byteNo+1])<<2 ] ;
    res[rc++] = '=';
  }
  
  res[rc]=0; // NULL TERMINATOR! ;)
  return res ;
}

Solver3d::Solver3d(){

    // lendo os parametros
    cout << "Lendo os parametros..." << endl;
    this->leParametros("parametros3d.txt");

    cout << "Parametros lidos com sucesso!" << endl;
    cout << endl;
    cout << "X  = " << this->d.X  << "m" << endl;
    cout << "Y  = " << this->d.Y  << "m" << endl;
    cout << "Z  = " << this->d.Z  << "m" << endl;
    cout << "T  = " << this->d.T  << "s" << endl;
    cout << "dx = " << this->d.dx << "m" << endl;
    cout << "dt = " << this->d.dt << "s" << endl;
    cout << "Nx = " << this->d.Nx << endl;
    cout << "Ny = " << this->d.Ny << endl;
    cout << "Nz = " << this->d.Nz << endl;
    cout << "Nt = " << this->d.Nt << endl;
    cout << "xs = " << this->d.xs << "m" << endl;
    cout << "ys = " << this->d.ys << "m" << endl;
    cout << "zs = " << this->d.zs << "m" << endl;
    cout << endl;

    cout << "Inicializando vetores..." << endl << endl;
    // inicializa vetores
    this->u_current = new Grid3d(this->d.Nx, this->d.Ny, this->d.Nz);
    this->u_next    = new Grid3d(this->d.Nx, this->d.Ny, this->d.Nz);

}

Solver3d::~Solver3d(){}

void Solver3d::leParametros(string nome){

    // * funcao que le os parâmetros do dominio 3D

    ifstream myfile;

    myfile.open("../" + nome);

    if(myfile.is_open()){

        // dimensoes do dominio
        myfile >> this->d.X;
        myfile >> this->d.Y;
        myfile >> this->d.Z;
        myfile >> this->d.T;

        // largura da malha
        myfile >> this->d.dx;
        this->d.dy = this->d.dx;
        this->d.dz = this->d.dx;
        myfile >> this->d.dt;

        // posicao da fonte
        myfile >> this->d.fcorte;
        myfile >> this->d.xs;
        myfile >> this->d.ys;
        myfile >> this->d.zs;

        // numero de iteracoes
        this->d.Nx = this->d.X/this->d.dx;
        this->d.Ny = this->d.Y/this->d.dy;
        this->d.Nz = this->d.Z/this->d.dz;
        this->d.Nt = this->d.T/this->d.dt;

        // velocidade da onda
        this->d.c = 2200;

        // this->d.vel = new Grid2d(this->d.Nz, this->d.Nx);
        // int v;

        // // matriz de velocidades
        // for (int j = 0; j < this->d.Nz; j++){
        //     for (int i = 0; i < this->d.Nx; i++){
        //         myfile >> v;
        //         this->d.vel->set(j, i, v);
        //     }
        // }

        myfile.close();

    } else {
        cerr << "Falha ao abrir arquivo de parametros" << endl;
        exit;
    }
    
}

void Solver3d::salvaVTI(Dominio d, Grid3d* u, string nomeDoArq, string info){


    // * Função que gera um arquivo vtk ImageData para o ParaView

    ofstream myfile;

    // cout << "Gerando arquivo data" << to_string(k/modk) << ".vti" << "..." << endl;

    myfile.open("../" + nomeDoArq + ".vti");

    if(myfile.is_open()){

        myfile << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
        myfile << "  <ImageData WholeExtent= \"" <<  STENCIL << " " << d.Nx - 1 - STENCIL << " " << STENCIL << " " << d.Ny - 1 - STENCIL << " " << STENCIL << " " << d.Nz - 1 - STENCIL << "\" ";
        myfile << "Origin = \"" << STENCIL << " " << STENCIL << " " << d.Nz - 1 - STENCIL << "\" ";
        myfile << "Spacing = \"" << d.dx << " " << d.dy << " " << d.dz << "\">\n";
        myfile << "    <Piece Extent = \"" << STENCIL << " " << d.Nx - 1 - STENCIL << " " << STENCIL << " " << d.Ny - 1 - STENCIL << " " << STENCIL << " " << d.Nz - 1 - STENCIL << "\">\n";
        myfile << "      <PointData Scalars=\"" + info + "\">\n";
        myfile << "        <DataArray type=\"Float32\" Name=\"" + info + "\" format=\"ascii\">\n";

        for (int k = STENCIL; k < d.Nz - STENCIL; k++){
            for (int j = STENCIL; j < d.Ny - STENCIL; j++){
                for (int i = STENCIL; i < d.Nx - STENCIL; i++){
                    myfile << u->get(k, j, i) << " ";
                }
            }
        }

        myfile << "\n        </DataArray>";
        myfile << "\n      </PointData>";
        myfile << "\n    </Piece>";
        myfile << "\n  </ImageData>";
        myfile << "\n</VTKFile>";

    } else {
        cout << "Erro na gravação do arquivo " << nomeDoArq << ".vti" << endl;
    }


}

void Solver3d::salvaVTIBinary(Dominio d, Grid3d* u, string nomeDoArq, string info){

    // * Função que gera um arquivo vtk ImageData para o ParaView

    ofstream myfile;

    myfile.open("../" + nomeDoArq + ".vti", ios::binary);

    if(myfile.is_open()){

        myfile << "<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\">\n";
        myfile << "  <ImageData WholeExtent= \"" <<  0 << " " << d.Nx-1 << " " << 0 << " " << d.Ny-1 << " " << 0 << " " << d.Nz-1 << "\" ";
        myfile << "Origin = \"" << 0 << " " << 0 << " " << d.Nz-1 << "\" ";
        myfile << "Spacing = \"" << d.dx << " " << d.dy << " " << d.dz << "\">\n";
        myfile << "    <Piece Extent = \"" << 0 << " " << d.Nx-1 << " " << 0 << " " << d.Ny-1 << " " << 0 << " " << d.Nz-1 << "\">\n";
        myfile << "      <PointData Scalars=\"" + info + "\">\n";
        myfile << "        <DataArray type=\"Float32\" Name=\"" + info + "\" format=\"appended\" offset=\"0\">\n";

        myfile << "        </DataArray>\n";
        myfile << "      </PointData>\n";
        myfile << "    </Piece>\n";
        myfile << "  </ImageData>\n";
        myfile << "  <AppendedData encoding=\"raw\">\n_";
        // int result_size;
        // char *encoding = base64((char *)u->firstptr(), sizeof(float) * u->getSize(), &result_size);
        // myfile << result_size;
        // myfile.write(encoding, result_size);
        int size = u->getSize()*sizeof(float);
        myfile << size;
        myfile.write((char *) u->firstptr(), size);
        myfile << "\n  </AppendedData>\n";
        myfile << "</VTKFile>\n";

    } else {
        cout << "Erro na gravação do arquivo " << nomeDoArq << ".vti" << endl;
    }


}


float Solver3d::pulso3d(int x, int y, int z, int t){

    // *funcao que simula um pulso sismico na posicao (xs, zs)

    if (x*this->d.dx != this->d.xs || y*this->d.dy != this->d.ys || z*this->d.dz != this->d.zs){
        return 0;
    } 

    float td = t*this->d.dt - ((2*sqrt(M_PI))/this->d.fcorte);  
    float fc = (this->d.fcorte/(3*sqrt(M_PI)));
    float r = (1.0 - 2.0 * M_PI * pow(M_PI * fc * td, 2))/pow(M_E, M_PI*pow((M_PI*fc*td), 2));
    // cout << r << endl;
    return r;

}

void Solver3d::mdf3d(Grid3d* u1, Grid3d* u2, int t){

    // * Calcula u2 a partir de u1 por DF

    float val;
    float courantNumber = d.dt*d.c/d.dx;
    float const1 = (pow(courantNumber, 2)/12.0);
    float const2 = pow(d.c*d.dt, 2);

    #pragma omp parallel for collapse(2) private(val)
    for (int k = STENCIL; k <= d.Nz - STENCIL; k++){
        for (int j = STENCIL; j <= d.Ny - STENCIL; j++){
            for (int i = STENCIL; i <= d.Nx - STENCIL; i++){

                val = 
                const1 *
                (
                -1*(u1->get(k,j,i-2) + u1->get(k,j,i+2) + u1->get(k,j-2,i) + u1->get(k,j+2,i) + u1->get(k-2,j,i) + u1->get(k+2,j,i)) +
                16*(u1->get(k,j,i-1) + u1->get(k,j,i+1) + u1->get(k,j-1,i) + u1->get(k,j+1,i) + u1->get(k-1,j,i) + u1->get(k+1,j,i)) -
                90*u1->get(k, j, i)
                ) + 
                2*u1->get(k, j, i) - u2->get(k, j, i) + this->pulso3d(i, j, k, t);

                u2->set(k, j, i, val);

            }
        }
    }

}


void Solver3d::mdf3d(Grid3d* u1, Grid3d* u2){

    // * Calcula u2 a partir de u1 por DF

    float val;
    float courantNumber = d.dt*d.c/d.dx;
    float const1 = (pow(courantNumber, 2)/12.0);

    for (int k = STENCIL; k < d.Nz - STENCIL; k++){
        for (int j = STENCIL; j < d.Ny - STENCIL; j++){
            for (int i = STENCIL; i < d.Nx - STENCIL; i++){

                val = 
                const1 *
                (
                -1*(u1->get(k,j,i-2) + u1->get(k,j,i+2) + u1->get(k,j-2,i) + u1->get(k,j+2,i) + u1->get(k-2,j,i) + u1->get(k+2,j,i)) +
                16*(u1->get(k,j,i-1) + u1->get(k,j,i+1) + u1->get(k,j-1,i) + u1->get(k,j+1,i) + u1->get(k-1,j,i) + u1->get(k+1,j,i)) -
                90*u1->get(k, j, i)
                ) + 
                2*u1->get(k, j, i) - u2->get(k, j, i);

                u2->set(k, j, i, val);

            }
        }
    }

}

void Solver3d::reynolds(){

    float courantNumber = d.dt * d.c/d.dx;

    // plano x = 0 e x = Nx
    #pragma omp parallel for collapse(2)
    for (int k = STENCIL; k <= d.Nz - STENCIL; k++){
        for(int j = STENCIL; j <= d.Ny - STENCIL; j++){
            u_next->set(k, j, STENCIL, u_current->get(k, j, STENCIL) + courantNumber*(u_current->get(k, j,STENCIL+1) - u_current->get(k, j, STENCIL)));
            u_next->set(k, j, d.Nx-STENCIL, u_current->get(k, j, d.Nx-STENCIL) - courantNumber*(u_current->get(k, j, d.Nx-STENCIL) - u_current->get(k, j, d.Nx-STENCIL-1)));
        }
    }

    // plano y = 0 e y = Ny
    #pragma omp parallel for collapse(2)
    for (int k = STENCIL; k <= d.Nz - STENCIL; k++){
        for(int i = STENCIL; i <= d.Nx - STENCIL; i++){
            u_next->set(k, STENCIL, i, u_current->get(k, STENCIL, i) + courantNumber*(u_current->get(k, STENCIL+1, i) - u_current->get(k, STENCIL, i)));
            u_next->set(k, d.Ny-STENCIL, i, u_current->get(k, d.Ny-STENCIL, i) - courantNumber*(u_current->get(k, d.Ny-STENCIL, i) - u_current->get(k, d.Ny-STENCIL-1, i)));
        }
    }

    // plano z = 0 e z = Nz
    #pragma omp parallel for collapse(2)
    for (int j = STENCIL; j <= d.Ny - STENCIL; j++){
        for(int i = STENCIL; i <= d.Nx - STENCIL; i++){
            u_next->set(STENCIL, j, i, u_current->get(STENCIL, j, i) + courantNumber*(u_current->get(STENCIL+1, j, i) - u_current->get(STENCIL, j, i)));
            u_next->set(d.Nz-STENCIL, j, i, u_current->get(d.Nz-STENCIL, j, i) - courantNumber*(u_current->get(d.Nz-STENCIL, j, i) - u_current->get(d.Nz-STENCIL-1, j, i)));
        }
    }

}

float Solver3d::atenuacao(float x, int borda){

    // * Função aplicada nas fronteiras para reduzir a amplitude da onda

    float fat = 0.0035;
    return exp(-(pow(fat*(borda - x), 2)));

}

void Solver3d::aplicaAmortecimento(){

    // * Percorre os planos das matrizes atual e proxima aplicando uma função de atenuação

    int borda = 20 + STENCIL;

    #pragma omp parallel 
    {

        // plano superior z = 0
        #pragma omp for collapse(2)
        for (int k = STENCIL; k <= borda; k++){
            for(int j = STENCIL; j <= d.Ny - STENCIL; j++){
                for(int i = STENCIL; i <= d.Nx - STENCIL; i++){

                    u_current->set(k, j, i, u_current->get(k, j, i)*atenuacao(k, borda));
                    u_next->set(k, j, i, u_next->get(k, j, i)*atenuacao(k, borda));
                    
                }
            }
        }

        // plano inferior z = Nz
        #pragma omp for collapse(2)
        for (int k = d.Nz - borda; k <= d.Nz - STENCIL; k++){
            for(int j = STENCIL; j <= d.Ny - STENCIL; j++){
                for(int i = STENCIL; i <= d.Nx - STENCIL; i++){

                    u_current->set(k, j, i, u_current->get(k, j, i)*atenuacao(d.Nz - k, borda));
                    u_next->set(k, j, i, u_next->get(k, j, i)*atenuacao(d.Nz - k, borda));
                    
                }
            }
        }

        // plano lateral y = 0
        #pragma omp for collapse(2)
        for (int k = STENCIL; k <= d.Nz - STENCIL; k++){
            for(int j = STENCIL; j <= borda; j++){
                for(int i = STENCIL; i <= d.Nx - STENCIL; i++){

                    u_current->set(k, j, i, u_current->get(k, j, i)*atenuacao(j, borda));
                    u_next->set(k, j, i, u_next->get(k, j, i)*atenuacao(j, borda));
                    
                }
            }
        }

        // plano lateral y = Ny
        #pragma omp for collapse(2)
        for (int k = STENCIL; k <= d.Nz - STENCIL; k++){
            for(int j = d.Ny - borda; j <= d.Ny - STENCIL; j++){
                for(int i = STENCIL; i <= d.Nx - STENCIL; i++){

                    u_current->set(k, j, i, u_current->get(k, j, i)*atenuacao(d.Ny - j, borda));
                    u_next->set(k, j, i, u_next->get(k, j, i)*atenuacao(d.Ny - j, borda));
                    
                }
            }
        }

        // plano lateral x = 0
        #pragma omp for collapse(2)
        for (int k = STENCIL; k <= d.Nz - STENCIL; k++){
            for(int j = STENCIL; j <= d.Ny - STENCIL; j++){
                for(int i = STENCIL; i <= borda; i++){

                    u_current->set(k, j, i, u_current->get(k, j, i)*atenuacao(i, borda));
                    u_next->set(k, j, i, u_next->get(k, j, i)*atenuacao(i, borda));
                    
                }
            }
        }

        // plano lateral x = Nx
        #pragma omp for collapse(2)
        for (int k = STENCIL; k <= d.Nz - STENCIL; k++){
            for(int j = STENCIL; j < d.Ny - STENCIL; j++){
                for(int i = d.Nx - borda; i <= d.Nx - STENCIL; i++){

                    u_current->set(k, j, i, u_current->get(k, j, i)*atenuacao(d.Nx - i, borda));
                    u_next->set(k, j, i, u_next->get(k, j, i)*atenuacao(d.Nx - i, borda));
                    
                }
            }
        }

    }

    

}

void Solver3d::solve(){

    int modk = 100;
    int t; // iterador temporal

    auto inicio = chrono::high_resolution_clock::now();

    // laço temporal com fonte
    cout << "Iniciando laço temporal..." << endl;
    for (t = 0; t < d.Nt; t += 2){
        
        // * calcula u_next
        this->mdf3d(u_current, u_next, t);
        // reynolds();
        aplicaAmortecimento();

        // * gera arquivo de dados a cada 100 iteracoes em k
        if (t % modk == 0){
            string nomeDoArq = "data" + to_string(t/modk);
            this->salvaVTIBinary(d, u_current, nomeDoArq, "P-Wave");
        }
        
        // * calcula u_current
        this->mdf3d(u_next, u_current, t + 1);
        aplicaAmortecimento();

    }

    cout << "\nArquivos gerados com sucesso!" << endl;

    auto final = chrono::high_resolution_clock::now();
    chrono::duration<double> intervalo = final - inicio;
    cout << "\nTempo decorrido: " << intervalo.count() << "s\n";

}