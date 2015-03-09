b=zeros(18,5,72);
for i=1:18
    for k=1:5
        sum0(i,k)=0;
        for j=1:72
            if S(i,j)==k
                sum0(i,k)=sum0(i,k)+1;
            end
        end
    end
end
for i=1:18
    j=1;
    for k=1:5
        s=sum0(i,k);
        while s~=0
            b(i,k,j)=sum0(i,k)/72;
            j=j+1;
            s=s-1;
        end
    end
end

for i=1:5
    b2(i,:)=b1(i,:)/aa(i);
end
B=b2;