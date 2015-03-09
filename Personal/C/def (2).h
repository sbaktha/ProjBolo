#define T 19
#define N 6
#define K 3 
          //for matrix B
double a[N][N]={{0.2,0.2,0.15,0.15,0.1,0.1},{0,0.2,0.1,0.25,0.25,0.1},{0,0,0.15,0.15,0.2,0.2},{0,0,0,0.25,0.3,0.45},{0,0,0,0,0.62,0.38},{0,0,0,0,0,1}};   // forward computation .
double b[N][K]={{0.4,0.4,0.2},{0.25,0.45,0.3},{0.2,0.35,0.45},{0.2,0.3,0.5},{0.6,0.2,0.2},{0.1,0.4,0.5}};
double Pi[N]={0.4,0.3,0.3}, sum;
int Ob[T]={1,2,1,2,1,0,1,1,1,0,2,1,0,0,1,2,2,1,0};   // possible values in the sequence .
double Beta[T][N]={0.0};   // Backward computation .
double Alpha[T][N]={0.0};
int i,j,t;//  T is for sequence length ;; N is for number of states 
double ZI[T][N][N]={0.0},nu=0.0;
double Gamma[T][N]={0.0};
double E_T[N]={0.0};
double E_I_J[N]={0.0};
double E_Pi[N]={0.0};
double E_A[N][N]={0.0};
double N_E_A[N][N]={0.0};
double E_B[N][N]={0.0};
double sum1[K]={0}; 
double p_v[N]={0.0};
int status[N]={0};
int m,n,tt;
void forward();
void backward();
void Bw_algo();
void P_V();
