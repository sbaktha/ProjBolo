Num1A=zeros(1,NoOfOb);
Den1A=zeros(1,NoOfOb);
Num1B=zeros(1,NoOfOb);
Den1B=zeros(1,NoOfOb);
NN_E_A=zeros(N,N);
New_A=zeros(N,N);
New_B=zeros(N,K);
NN_E_B=zeros(N,K);

for i=1:N
    for j=1:N
        for kk=1:NoOfOb
            Ob=Observtot(order(kk),:);
            sum0=0;
            sum1=0;
            for t=1:T-1
                sum0 = sum0 + (Alphastore(kk,t,i) *astore(kk,i,j) *bstore(kk,j,Ob(t+1)) *Betastore(kk,t+1,j));
                sum1=sum1+Alphastore(kk,t,i)*Betastore(kk,t,i);
            end
            Num1A(kk)=(sum0)/FinalProb(kk);
            Den1A(kk)=(sum1)/FinalProb(kk);
        end
        New_A(i,j)=sum(Num1A(:))/sum(Den1A(:));
    end
end

for i=1:N
    for l=1:K
        for kk=1:NoOfOb
            Ob=Observtot(order(kk),:);
            sum0=0;
            sum1=0;
            for t=1:T-1
                if Ob(t)==l
                    sum0=sum0+Alphastore(kk,t,i)*Betastore(kk,t,i);
                end
                sum1=sum1+Alphastore(kk,t,i)*Betastore(kk,t,i);
            end
            Num1B(kk)=sum0/FinalProb(kk);
            Den1B(kk)=sum1/FinalProb(kk);
        end
        New_B(i,l)=sum(Num1B(:))/sum(Den1B(:));
    end
end
% counterb=zeros(1,N);
% for i=1:N
%     for j=1:K
%         if New_B(i,j)==0
%             counterb(i)=counterb(i)+1;
%             New_B(i,j)=0.0001;
%         end
%     end
% end
% for i=1:N
%     New_B(i,:)=New_B(i,:)/(sum(New_B(i,:))+counterb(i)*0.0001);
% end

