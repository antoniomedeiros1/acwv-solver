class Matriz3D{

    public:

        int l; //largura
        int a; //altura
        int p; //profundidade
        float* mat; //a matriz será representada por um vetor de tamanho l*a*p

        Matriz3D(int x, int z, int t){

            l = x;
            a = z;
            p = t;
            mat = new float[l*a*p];

            for (int i = 0; i < l*a*p; i++){
                mat[i] = 0;
            }

        }

        ~Matriz3D(){
            delete mat;
        }
        
        float get(int i, int j, int k){
            if (i < 0 || i > l || j < 0 || j > a || k < 0 || k > p){
                return 0;
            }
            return mat[k*(l*a) + j*l + i];
        }

        void set(int i, int j, int k, float val){
            if (!(i < 0 || i > l || j < 0 || j > a || k < 0 || k > p)){
                mat[k*(l*a) + j*l + i] = val;
            }
        }
        

};