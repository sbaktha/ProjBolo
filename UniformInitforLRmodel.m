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

        
                
