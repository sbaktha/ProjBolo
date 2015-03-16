Num1A=zeros(N,N);
Den1A=zeros(N,N);
Num1B=zeros(N,N);
Den1B=zeros(N,N);
NN_E_A=zeros(N,N);
New_A=zeros(N,N);
New_B=zeros(N,K);
NN_E_B=zeros(N,K);

for i=1:N
    for j=1:N
        for kk=1:NoOfOb
            sum0=zeros(1,T);
            sum1=zeros(1,T);
            for t=1:T-1
                sum0(t) = sum0(t) + (Alpha(t,i) *a(i,j) *b(j,Ob(t+1)) *Beta(t+1,j));
                sum1(t)=sum1(t)+Alpha(t,i)*Beta(t,i);
            end
            Num1A(kk)=Num1A(kk)+sum(sum0(:));
            Den1A(kk)=Den1A(kk)+sum(sum1(:));
            
        end
    end
end



if Ob(t)==L
    sum2(t)=sum2(t)+Alpha(t,i)*Beta(t,i);
end