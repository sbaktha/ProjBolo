% Definitions
% global T;global N;global K;global a;global b;global Pi;global Ob;global Beta;global Alpha;
% global i;global j;global t;global ZI;global nu;global Gamma;global E_T;global E_I_J;global E_Pi;global N_E_A;global E_A;global E_B;global sum1;global sum;
% global p_v;global m;global n;global status;global tt;
T = 72;
N = 5;
K = 72;
A  = [0.2,0.2,0.15,0.15,0.1,0.1 ; 0,0.2,0.1,0.25,0.25,0.1 ; 0,0,0.15,0.15,0.2,0.2 ; 0,0,0,0.25,0.3,0.45 ; 0,0,0,0,0.62,0.38 ; 0,0,0,0,0,1];
B  = [0.4,0.4,0.2 ; 0.25,0.45,0.3 ; 0.2,0.35,0.45 ; 0.2,0.3,0.5 ; 0.6,0.2,0.2 ; 0.1,0.4,0.5];
Pi = [1,0,0,0,0,0];
sum0=0;
%Ob=randi(30,1,72);
Ob= sort(randi(72,1,72));
Beta=zeros(T,N);
Alpha=zeros(T,N);
i=0;
j=0;
t=0;
ZI=zeros(T,N,N);
nu=0.0;
Gamma=zeros(T,N);
E_T=zeros(1,N);
E_I_J=zeros(1,N);
E_Pi=zeros(1,N);
E_A=zeros(N,N);
N_E_A=zeros(N,N);
E_B=zeros(N,K);
sum1=zeros(K);
p_v=zeros(N);
status=zeros(1,N);
m=0;
n=0;
tt=0;

%Forward Algorithm
for i=1:N
    Alpha(1,i)=Pi(i) * B(i,Ob(1));
end

for t=1:T-1
    for j=1:N
        sum0=0;
        for i=1:N
            sum0= sum0+ Alpha(t,i)*A(i,j);
        end
        Alpha(t+1,j)=sum0 * B(j,Ob(t+1));
    end
end
disp('The forward matrix is:');

for i=1:T
    for j=1:N
        fprintf('%.8f',Alpha(i,j));
        fprintf('    ');
    end
    fprintf('\n');
end
fprintf('\n');

%Backward Algorithm
for i=1:N
    Beta(T,i)=1;
end

for t=T-1:-1:1
    for i=1:N
        sum0=0;
        for j=1:N
            sum0=sum0+(A(i,j)*Beta(t+1,j)*B(j,Ob(t+1)));
        end
        Beta(t,i)=sum0;
    end
end
disp('The backward matrix is:');

for i=1:T
    for j=1:N
        fprintf('%.8f',Beta(i,j));
        fprintf('    ');
    end
    fprintf('\n');
end
fprintf('\n');

%Baum-Welch Algorithm
kk=0;
sum2=0;
%Calculation of ZI values
for t=1:T-1
    for i=1:N
        for j=1:N
            nu=Alpha(t,i)*B(j,Ob(t+1))*Beta(t+1,j)*A(i,j);
            sum0=0;
            for m=1:N
                for n=1:N
                    sum0 = sum0 + (Alpha(t,m) *a(m,n) *B(n,Ob(t+1)) *Beta(t+1,n));
                end
            end
            ZI(t,i,j) = nu/sum0;
        end
    end
end
% disp('The ZI matrix is:');
% disp(ZI);

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
disp('The Gamma matrix is:');               
for i=1:T
    for j=1:N
        fprintf('%.8f',Gamma(i,j));
        fprintf('    ');
    end
    fprintf('\n');
end
fprintf('\n');

%Expected number of transistions from state i
for i=1:N
    sum0=0;
    for t=1:T-1
        sum0= sum0 + Gamma(t,i);
    end
    E_T(i)=sum0;
end
disp('Expected no of transitions from the states:');
for i=1:N
fprintf('%.4f',E_T(i));
fprintf('\n');
end

%Expected number of transitions from node i to node j
for i=1:N
    for j=1:N
        sum0=0;
        for t=1:T-1
            sum0= sum0+ZI(t,i,j);
        end
        E_I_J(i)=sum0;                   % may be a mistake (already mentioned in the C version)..............
        fprintf('Expected no of transitions from the state %d to state %d:', i,j);
        fprintf('%.4f \n',E_T(i));
    end
end

%Computing estimated values for Pi ,A and B.
for i=1:N
    E_Pi(i)= Gamma(1,i);    % E_Pi(i)= Gamma(0,(i));
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
    for t=1:T  %to traverse the observation sequence...
        for kk=1:K
            if(Ob(t) == kk)  % here one for loop will come
                sum2 = sum2+ Gamma(t,j); % overall sum ..........
                sum1(kk)= sum1(kk) + Gamma(t,j);
                break;
            end
        end
    end
    for kk=1:K
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
disp(N);
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
end

tt=1;
for i=1:N
    if(p_v(i)*100 >= 40.0)
        status(tt)=i;
        tt= tt +1;
    else
        status(i)=0;
    end
end
disp('The status during the transition is:');
disp(status);
fprintf('\n');



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
                sum3= sum3 + A(i,j);
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
                N_E_A(i,j)=A(i,j);
            end
        end
    end
end
disp('After Normalization:');
fprintf('\n');
disp('The estimated state transition matrix is:');
for i=1:N
    for j=1:N
        fprintf('%.8f',E_A(i,j));
        fprintf('    ');
    end
    fprintf('\n');
end
fprintf('\n');
disp('The estimated probability matrix is:');
for i=1:N
    for j=1:K
        fprintf('%.8f',E_B(i,j));
        fprintf('    ');
    end
    fprintf('\n');
end

