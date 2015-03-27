
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
                tofindmax(i)=Delta(t-1,i)*a(i,j);
            end
            [Delta(t,j),Shi(t,j)]=max(tofindmax);
            clear tofindmax;
        end
    end
    for i=1:N
        tofindmax(i)=Delta(T,i);
    end
    [Pnew,StatenewViter(T)]=max(tofindmax);
    
    for t=T-1:-1:1
        StatenewViter(t)=Shi(t+1)*StatenewViter(t+1);
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
    b(E_B==0)=0.0001;    
    for i=1:4
        b(i,:)=bsxfun(@rdivide,b(i,:),sum(b(i,:)));
    end
end
%mlop(xxx,:);
