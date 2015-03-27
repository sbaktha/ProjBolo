def;
global Alpha;global Beta;global ZI;global Gamma; global E_A;global E_Pi;global E_B;global status;global N_E_A;
Alpha=forw();
Beta=backw();
[ZI,Gamma,E_Pi,E_A,E_B]=Bw_algo();
[status]=P_V();
[N_E_A, E_A]= normal();
