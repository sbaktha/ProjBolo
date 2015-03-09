
N=st;
for i=1:N
    pii(i)=0;
    pc(i)=0;
    for z=1:ex
        if(S(z,1)==i)
            pii(i)=pii(i)+1;
        end
    end
    Pi(i)=pii(i)/ex;
end

A=zeros(N);
AA=zeros(N);
for i=1:N
    for j=1:N
        AA(i,j)=0;
        for y=1:ex
            for x=1:wind-1
                if((S(y,x)==i)&&(S(y,(x+1))==j))
                    AA(i,j)=AA(i,j)+1;
                end
            end
        end
        te=sum(AA,2);
    end
    A(i,:)=AA(i,:)/te(i);
end

for k=1:st
    for j=1:symb
        sum0(k,j)=0;
        for i=1:ex
            for m=1:72
                if (S(i,m)==k && Observ(i,m)==j)
                    sum0(k,j)=sum0(k,j)+1;
                end
            end
        end
    end
end

b2=zeros(size(sum0));
sumstate1=sum(sumstate,2);
for i=1:st
    b2(i,:)=sum0(i,:)/sumstate1(i);
end

%%%%------------Viterbi for state optimizxation ----------%%%%
kluster=wind;
clear b;
b=zeros(ex,st,kluster);
for i=1:ex
    for k=1:st
        summ(i,k)=0;
        for j=1:kluster
            if S(i,j)==k
                summ(i,k)=summ(i,k)+1;
            end
        end
    end
end
for i=1:ex
    j=1;
    for k=1:st
        s=summ(i,k);
        while s~=0
            b(i,k,j)=summ(i,k)/kluster;
            j=j+1;
            s=s-1;
        end
    end
end
