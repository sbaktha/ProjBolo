%Forward Algorithm
for i=1:N
    Alpha(1,i)=Pi(i) * B(i,Ob(1));
end
for t=1:T-1
    for j=1:N
        sum=0;
        for i=1:N
            sum= sum+ Alpha(t,i)*A(i,j);
        end
        Alpha(t+1,j)=sum * B(j,Ob(t+1));
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
        sum=0;
        for j=1:N
            sum=sum+(A(i,j)*Beta(t+1,j)*B(j,Ob(t+1)));
        end
        Beta(t,i)=sum;
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
            nu=Alpha(t,i)*A(i,j)*B(j,Ob(t+1))*Beta(t+1,j);
            sum=0;
            for z=1:N
                for n=1:N
                    sum = sum + (Alpha(t,z) *A(z,n) *B(n,Ob(t+1)) *Beta(t+1,n));
                end
            end
            ZI(t,i,j) = nu/sum;
        end
    end
end
% disp('The ZI matrix is:');
% disp(ZI);

%Gamma computation
for t=1:T
    for i=1:N
        sum=0;
        for j=1:N
            sum = sum+ZI(t,i,j);
        end
        Gamma(t,i)=sum;
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
    sum=0;
    for t=1:T-1
        sum= sum + Gamma(t,i);
    end
    E_T(i)=sum;
end
disp('Expected no of transitions from the states:');
for i=1:N
fprintf('%.4f',E_T(i));
fprintf('\n');
end

%Expected number of transitions from node i to node j
for i=1:N
    for j=1:N
        sum=0;
        for t=1:T-1
            sum= sum+ZI(t,i,j);
        end
        E_I_J(i)=sum;                   % may be a mistake (already mentioned in the C version)..............
        fprintf('Expected no of transitions from the state %d to state %d:', i,j);
        fprintf('%.4f \n',E_I_J(i));
    end
end

%Computing estimated values for Pi ,A and B.
for i=1:N
    E_Pi(i)= Gamma(1,i);    % E_Pi(i)= Gamma(0,(i));
end

for i=1:N
    for j=1:N
        sum=0;
        nu=0;
        for t=1:T-1
            sum=sum+ZI(t,i,j);
            nu=nu+Gamma(t,i);
        end
        E_A(i,j) = (sum / nu)  ;
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




W_S_J_E_A=E_A;

%Computing the matrix B
for j=1:N     % number of states
    sum2=0;
    for kk=1:K
        sum1(kk)=0;
    end
    for t=1:T  %to traverse the observation sequence...
        for kk=1:K
            if(Ob(t) == kk)  % here one for loop will come
                sum2 = sum2+ Gamma(t,j);
                % overall sum ..........
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

W_S_J_E_B=E_B;

%probability of visit
sum = 0;
disp(' The probability of the node being visited during the training phase');
p_v=zeros(5,1);
for i=1:N
    if(i==1)
        p_v(i)=E_Pi(i);
    else
        sum=0;
        for j=1:(i-1)
            sum= sum + p_v(j)*(E_A(j,i)/(1-E_A(j,j)) );
        end
    end
    p_v(i)= sum + (E_Pi(i));
    disp(i);
    disp(p_v(i));
end

tt=1;u=1;
for i=1:N
        if(p_v(i)*100 >= 70.0)
            status(tt)=i;
            tt= tt +1;
        else
            skippedNode(1,u)=i;
            u=u+1;
            status(tt)=0;
            tt=tt+1;
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
        fprintf('%.8f',N_E_A(i,j));
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



diffEA=W_S_J_E_A-N_E_A;
% diffPi=W_S_J_E_Pi-N_E_Pi;
