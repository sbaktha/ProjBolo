def;
global Alpha;global Beta;global ZI;global Gamma; global E_A;global E_Pi;global E_B;global status;global N_E_A;
Alpha=forw();
Beta=backw();
[u,v,w,x,y]=Bw_algo();
ZI=u;
Gamma=v;
E_Pi=w;
E_A=x;
E_B=y;
[sta]=P_V();
status=sta;
[nonerg, erg]= noramal();
N_E_A=nonerg;
E_A=erg;