% All the declared variables and arrays have been assigned 0 unless and
% otherwise stated as in the case of arrays 'a' and 'b'. The function
% declaration has been omitted as there is no declaration part in MATLAB.
global T;global N;global K;global a;global b;global Pi;global Ob;global Beta;global Alpha;
global i;global j;global t;global ZI;global nu;global Gamma;global E_T;global E_I_J;global E_Pi;global N_E_A;global E_A;global E_B;global sum1;global sum;
global p_v;global m;global n;global status;global tt;
T = 19;  % T=20;
N = 6;   %N=5;
K = 3;   %k=4;
a = [0.2,0.2,0.15,0.15,0.1,0.1 ; 0,0.2,0.1,0.25,0.25,0.1 ; 0,0,0.15,0.15,0.2,0.2 ; 0,0,0,0.25,0.3,0.45 ; 0,0,0,0,0.62,0.38 ; 0,0,0,0,0,1];
b = [0.4,0.4,0.2 ; 0.25,0.45,0.3 ; 0.2,0.35,0.45 ; 0.2,0.3,0.5 ; 0.6,0.2,0.2 ; 0.1,0.4,0.5];
Pi=[0.4,0.3,0.3,0,0,0];
sum=0;
Ob=[2,3,2,3,2,1,2,2,2,1,3,2,1,1,2,3,3,2,1];
Beta=zeros(T,N);
Alpha=zeros(T,N);
i=0;
j=0;
t=0;
ZI=zeros(T,N,N);
nu=0.0;
Gamma=zeros(T,N);
E_T=zeros(N);
E_I_J=zeros(N);
E_Pi=zeros(N);
E_A=zeros(N,N);
N_E_A=zeros(N,N);
E_B=zeros(N,N);
sum1=zeros(K);
p_v=zeros(N);
status=zeros(N);
m=0;
n=0;
tt=0;