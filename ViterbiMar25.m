New1StartingCode
UniformInitforLRmodel
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
NoOfOb=18;
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
    %SingleLoopCalcHMM;
    %SingleLoopMar15DiffAB;
    
    looptimes=25;
    for yyy=1:looptimes
        %b
        mlop(xxx,yyy)=0;
        for i=1:N
            Alpha(1,i)=Pi(i) * b(i,Ob(1));
        end
        
        for t=1:T-1
            for j=1:N
                sum0=0;
                for i=1:N
                    sum0= sum0+ Alpha(t,i)*a(i,j);
                end
                Alpha(t+1,j)=sum0 * b(j,Ob(t+1));
            end
        end
        %Alpha;
        ProbOgivenLfwd(yyy)=0;
        for i=1:N
            ProbOgivenLfwd(yyy)=ProbOgivenLfwd(yyy)+Alpha(T,i);
        end
        
        for i=1:N
            Beta(T,i)=1;
        end
        
        for t=T-1:-1:1
            for i=1:N
                sum0=0;
                for j=1:N
                    sum0=sum0+(a(i,j)*Beta(t+1,j)*b(j,Ob(t+1)));
                end
                Beta(t,i)=sum0;
            end
        end
        %Beta;
        
        ProbOgivenLbkwd(yyy)=0;
        for i=1:N
            ProbOgivenLbkwd(yyy)=ProbOgivenLbkwd(yyy)+Pi(i)*b(i,Ob(1))*Beta(1,i);
        end
        
        
        ZI=zeros(T,N,N);
        kk=0;
        %Calculation of ZI values
        sum0=zeros(1,T);
        for t=1:T-1
            for i=1:N
                for j=1:N
                    sum0(t) = sum0(t) + (Alpha(t,i) *a(i,j) *b(j,Ob(t+1)) *Beta(t+1,j));
                end
            end
        end
        
        for t=1:T-1
            for i=1:N
                for j=1:N
                    nu=Alpha(t,i)*b(j,Ob(t+1))*Beta(t+1,j)*a(i,j);
                    ZI(t,i,j) = nu/sum0(t);
                end
            end
        end
        % Gamma alternate
        %     sum00=zeros(1,T);
        %     for t=1:T
        %         for i=1:N
        %             sum00(t)=sum00(t)+Alpha(t,i)*Beta(t,i);
        %         end
        %     end
        for t=1:T
            sum00=0;
            for i=1:N
                sum00=sum00+Alpha(t,i)*Beta(t,i);
            end
            for i=1:N
                Gammanew(t,i)=Alpha(t,i)*Beta(t,i)/sum00;
            end
        end
        
        
        for i=1:N
            Delta(1,i)=Pi(i)*b(i,Ob(1));
            Shi(1,i)=0;
        end
        for t=2:T
            for j=1:N
                for i=1:N
                    tofindmaxshi(i)=Delta(t-1,i)*a(i,j);
                    tofindmaxdel(i)=Delta(t-1,i)*a(i,j);
                end
                Delta(t,j)=max(tofindmaxdel)*b(j,Ob(t));
                [temp, Shi(t,j)]=max(tofindmaxshi);
                %                 Delta(t,j)=Delta(t,j)*b(j,Ob(t));
                clear tofindmax;
            end
        end
        for i=1:N
            tofindmax(i)=Delta(T,i);
        end
        [Pnew,StatenewViter(T)]=max(tofindmax);
        
        for t=T-1:-1:1
            StatenewViter(t)=Shi(t+1,StatenewViter(t+1));
        end
        
        for t=1:T-1
            for i=1:N
                sum0=0;
                for j=1:N
                    sum0 = sum0+ZI(t,i,j);
                end
                Gamma(t,i)=sum0;
            end
        end
        
        
        Gamma(T,:)=Gammanew(T,:);
        
        for i=1:N
            sum0=0;
            for t=1:T-1
                sum0= sum0 + Gamma(t,i);
            end
            E_T(i)=sum0;
        end
        
        for i=1:N
            for j=1:N
                sum0=0;
                for t=1:T-1
                    sum0= sum0+ZI(t,i,j);
                end
                E_I_J(i,j)=sum0;
            end
        end
        
        for i=1:N
            E_Pi(i)= Gamma(1,i);
        end
        for i=1:N
            for j=1:N
                sum0=0;
                nu=0;
                for t=1:T-1
                    sum0=sum0+ZI(t,i,j);
                    nu=nu+Gamma(t,i);
                end
                E_A(i,j) = (sum0 / nu)  ;
            end
        end
        for j=1:N     % number of states
            sum2=0;
            for kk=1:K
                sum1(kk)=0;
            end
            for t=1:T  %to traverse the ob  servation sequence...
                sum2 = sum2+ Gamma(t,j);
                for kk=1:K
                    if(Ob(t) == kk)
                        sum1(kk)= sum1(kk) + Gamma(t,j);
                        break;
                    end
                end
            end
            for kk=1:K
                E_B(j,kk) = (sum1(kk))/sum2;
            end
        end
        
        mlop(xxx,yyy)=mlop(xxx,yyy)+ProbOgivenLbkwd(yyy);
        
        a=E_A;
        b=E_B;
        b(E_B<0.0001)=0.0001;
        for i=1:4
            b(i,:)=bsxfun(@rdivide,b(i,:),sum(b(i,:)));
        end
    end
    %mlop(xxx,:);
    
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
    N_E_B(N_E_B<0.0001)=[0.0001];
    for i=1:4
        N_E_B(i,:)=bsxfun(@rdivide,N_E_B(i,:),sum(N_E_B(i,:)));
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
    anstore(xxx,:,:)=N_E_A;
    bnstore(xxx,:,:)=N_E_B;
end
mlogprob=zeros(1,looptimes);
for j=1:looptimes
    for i=1:NoOfOb
        mlogprob(j)=mlogprob(j)+log(mlop(i,j))/ex;
    end
end
AijBjlDiffcalcMar23;


