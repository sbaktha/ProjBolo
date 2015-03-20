clear b;
b=b3;
a=A;
N = st;
K = symb;
sum0=0;
Ob=zeros(size(Observ(1,:)));
T=length(Ob);
N_E_A=zeros(N,N);
E_B=zeros(N,K);
N_E_B=zeros(N,K);
sum1=zeros(K);
p_v=zeros(N);
status=zeros(1,N);
nu=0.0;
obiter=size(Ob,1);
%order=randperm(18,18);
order=[1:18];
NoOfOb=1;
looptimes=1;
for xxx=1:NoOfOb
    Beta=zeros(T,N);
    Alpha=zeros(T,N);
    ZI=zeros(T,N,N);
    Gamma=zeros(T,N);
    E_T=zeros(1,N);
    E_I_J=zeros(N,N);
    E_Pi=zeros(1,N);
    E_A=zeros(N,N);
    E_B=zeros(N,K);
    Ob=Observ(order(xxx),:);
    SingleLoopCalcHMM;
    %SingleLoopMar15DiffAB;
    Alphastore(xxx,:,:)=Alpha;
    Betastore(xxx,:,:)=Beta;
    FinalProb(xxx)=ProbOgivenLfwd(end);
    %mlogprob(xxx)=sum(log(mlop(xxx,:))/obiter);
    %Normalization
    sum2 = 0;
    sum3 = 0;
    pp = 0;
    pp1 = 0;
    for i=1:N
        if(i==status(pp+1)) % status(pp)
            pp=pp+1;
            for j=1:N
                N_E_A(i,j)=E_A(i,j);
            end
        else
            sum3=0;
            sum2=0;
            pp1=0;
            for j=1:N
                if(j==status(pp1+1)) %status(pp1)
                    pp1=pp1+1;
                    sum3= sum3 + a(i,j);
                else
                    sum2= sum2 + E_A(i,j);
                end
            end
            
            pp1=0;
            for j=1:N
                if(j~=status(pp1+1))  %status(pp1)
                    N_E_A(i,j)=(1-sum3)*(E_A(i,j)/sum2);
                else
                    pp1=pp1+1;
                    N_E_A(i,j)=a(i,j);
                end
            end
        end
    end
    
    
    N_E_B=E_B;
    N_E_B(N_E_B==0)=[0.00001];
    if sum(N_E_B==0)~=0
        for i=1:N
            coun(i)=sum(N_E_B(i,:)==0.00001);
        end
        for i=1:N
            N_E_B(i,:)=N_E_B(i,:)/(sum(E_B(i,:))+coun(i));
        end
    end
    a=N_E_A;
    b=N_E_B;
    Pi=E_Pi;
    
    disp('After Normalization:');
    fprintf('\n');
    disp('The estimated state transition matrix is:');
    for i=1:N
        for j=1:N
            fprintf('%.8f',N_E_A(i,j));
            fprintf('    ');
        end
        fprintf('\n');
    end
    fprintf('\n');
    disp('The estimated probability matrix is:');
    for i=1:N
        for j=1:K
            fprintf('%.8f',N_E_B(i,j));
            fprintf('    ');
        end
        fprintf('\n');
    end
    astore(xxx,:,:)=E_A;
    bstore(xxx,:,:)=E_B;
end
mlogprob=zeros(1,looptimes);
for j=1:looptimes
    for i=1:NoOfOb
        mlogprob(j)=mlogprob(j)+log(mlop(i,j))/ex;
    end
end
AijBjlDiffcalcMar16;
