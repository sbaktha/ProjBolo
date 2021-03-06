clear b;
b=b2;
a=A;
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

for xxx=1:1
    mlogprob(xxx)=0;
    for yyy=1:obiter
        mlop(xxx,yyy)=0;
        if mod(yyy,18)~=0
            Ob=Observ(mod(yyy,18),:);
        else
            Ob=Observ(mod(yyy,18)+1,:);
        end
        %Fowward Algorithm
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
        %     disp('The forward matrix is:');
        %
        %     for i=1:T
        %         for j=1:N
        %             fprintf('%.8f',Alpha(i,j));
        %             fprintf('    ');
        %         end
        %         fprintf('\n');
        %     end
        %     fprintf('\n');
        ProbOgivenLfwd(yyy)=0;
        for i=1:N
            ProbOgivenLfwd(yyy)=ProbOgivenLfwd(yyy)+Alpha(T,i);
        end
        
        %Backward Algorithm
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
        %     disp('The backward matrix is:');
        %
        %     for i=1:T
        %         for j=1:N
        %             fprintf('%.8f',Beta(i,j));
        %             fprintf('    ');
        %         end
        %         fprintf('\n');
        %     end
        %     fprintf('\n');
        %
        ProbOgivenLbkwd(yyy)=0;
        for i=1:N
            ProbOgivenLbkwd(yyy)=ProbOgivenLbkwd(yyy)+Pi(i)*b(i,Ob(1))*Beta(1,i);
        end
        
        
        %     %Viterbi for state optimization
        %     for i=1:N
        %         delta(1,i)=Pi(i)*b(i,Ob(1));
        %         shi(1,i)=0;
        %     end
        %     for j=1:N
        %         for t=2:T
        %             delta(t,j)=(max(i,delta(t-1,i)*a(i,j)))*b(j,Ob(t));
        %             shi(t,j)=argmax(i,delta(t-1,i)*a(i,j)));
        %         end
        %     end
        %
        %     ProbOandOS=max(i,delta(T,i));
        %     xnew(T)=argmax(i,delta(T,i));
        %
        %     for t=T-1:-1:1
        %         xnew(t)=shi(t+1,xnew(t+1));
        %     end
        
        %%%%%%%%%%%%%%%%%%%
        
        %Baum-Welch Algorithm
        kk=0;
        %Calculation of ZI values
        for t=1:T-1
            for i=1:N
                for j=1:N
                    nu=Alpha(t,i)*b(j,Ob(t+1))*Beta(t+1,j)*a(i,j);
                    sum0=0;
                    for m=1:N
                        for n=1:N
                            sum0 = sum0 + (Alpha(t,m) *a(m,n) *b(n,Ob(t+1)) *Beta(t+1,n));
                        end
                    end
                    ZI(t,i,j) = nu/sum0;
                end
            end
        end
        
        %Gamma computation
        for t=1:T
            for i=1:N
                sum0=0;
                for j=1:N
                    sum0 = sum0+ZI(t,i,j);
                end
                Gamma(t,i)=sum0;
            end
        end
        %     disp('The Gamma matrix is:');
        %     for i=1:T
        %         for j=1:N
        %             fprintf('%.8f',Gamma(i,j));
        %             fprintf('    ');
        %         end
        %         fprintf('\n');
        %     end
        %     fprintf('\n');
        
        %Expected number of transistions from state i
        for i=1:N
            sum0=0;
            for t=1:T-1
                sum0= sum0 + Gamma(t,i);
            end
            E_T(i)=sum0;
        end
        %     disp('Expected no of transitions from the states:');
        %     for i=1:N
        %         fprintf('%.4f',E_T(i));
        %         fprintf('\n');
        %     end
        
        %Expected number of transitions from node i to node j
        for i=1:N
            for j=1:N
                sum0=0;
                for t=1:T-1
                    sum0= sum0+ZI(t,i,j);
                end
                E_I_J(i)=sum0;
                %             fprintf('Expected no of transitions from the state %d to state %d:', i,j);
                %             fprintf('%.4f \n',E_I_J(i));
            end
        end
        
        %Computing estimated values for Pi ,A and B.
        for i=1:N
            E_Pi(i)= Gamma(1,i);
        end
        for i=1:N
            for j=1:N
                sum0=0;
                nu=0;
                for t=1:T-1
                    sum0=sum0+ZI(t+1,i,j);
                    nu=nu+Gamma(t,i);
                end
                E_A(i,j) = (sum0 / nu)  ;
            end
        end
        disp('The estimated state transition matrix is:');
        for i=1:N
            for j=1:N
                fprintf('%.8f',E_A(i,j));
                fprintf('    ');
            end
            fprintf('\n');
        end
        fprintf('\n');
        
        %Computing the matrix B
        for j=1:N     % number of states
            sum2=0;
            for kk=1:K
                sum1(kk)=0;
            end
            for t=1:T  %to traverse the ob  servation sequence...
                for kk=2:K
                    if(Ob(t) == kk)
                        sum2 = sum2+ Gamma(t,j); % overall sum ..........
                        sum1(kk)= sum1(kk) + Gamma(t,j);
                        break;
                    end
                end
            end
            for kk=2:K
                E_B(j,kk) = (sum1(kk))/sum2;
            end
        end
        disp('The estimated probability matrix is:');
        for i=1:N
            for j=1:K
                fprintf('%.8f',E_B(i,j));
                fprintf('    ');
            end
            fprintf('\n');
        end
        fprintf('\n');
        
        %probability of visit
        sum0 = 0;
        disp(' The probability of the node being visited during the training phase');
        for i=1:N
            if(i==1)
                p_v(i)=E_Pi(i);
            else
                sum0=0;
                for j=1:(i-1)
                    sum0= sum0 + p_v(j)*(E_A(j,i)/(1-E_A(j,j)) );
                end
            end
            p_v(i)= sum0 + (E_Pi(i));
            disp(i);
            disp(p_v(i));
        end
        
        tt=1;
        for i=1:N
            if(p_v(i)*100 >= 35.0)
                status(tt)=i;
                tt= tt +1;
            else
                status(i)=0;
            end
        end
        disp('The status during the transition is:');
        disp(status);
        fprintf('\n');
        
        %%%%%%
        mlop(xxx,yyy)=mlop(xxx,yyy)+ProbOgivenLbkwd(yyy);
    end
    
    mlogprob(xxx)=sum(log(mlop(xxx,:))/obiter);
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
end