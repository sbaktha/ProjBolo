%Uniform Initialization
N=st;
K=symb;
Pi=zeros(1,N);
A=zeros(N);
B=zeros(N,K);

Pi(1)=1;
for i=2:N
    Pi(i)=0;
end

for i=1:N
    for j=1:N
        if j==i || j==i+1
                A(i,j)=0.5;
        end
    end
end

for i=1:N
    for j=1:K
        B(i,j)=1/K;
    end
end

a=A;
b=B;
N = st;
K = symb;
sum0=0;
Ob=zeros(size(Observ(1,:)));
T=length(Ob);
Beta=zeros(T,N);
Alpha=zeros(T,N);
ZI=zeros(T,N,N);
Gamma=zeros(T,N);
E_T=zeros(1,N);
E_I_J=zeros(1,N);
E_Pi=zeros(1,N);
E_A=zeros(N,N);
N_E_A=zeros(N,N);
E_B=zeros(N,K);
N_E_B=zeros(N,K);
sum1=zeros(K);
p_v=zeros(N);
status=zeros(1,N);
nu=0.0;
obiter=size(Ob,1);
order=randperm(18,10);

