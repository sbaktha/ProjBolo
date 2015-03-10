
temp=kluster-1;
for i=1:temp
    for j=1:st
        for k=1:st
            temp1=0;temp2=0;
            for l=1:ex
                if S(l,i)==j
                    temp2=temp2+1;
                    if S(l,i+1)==k
                        temp1=temp1+1;
                    end
                end
            end
            if temp2==0
                aaa(i,j,k)=0;
            else
                aaa(i,j,k)=temp1/temp2;
            end
        end
    end
end
a1=aaa;

for i=1:st
    temp=0;
    for j=1:ex
        if S(j,1)==i
            temp=temp+1;
        end
    end
    py(i)=(temp/ex);
end
%%%%%%%%


for l=1:ex
    b1(:,:)=b(l,:,:);
    for i=1:st
        del(1,i)=py(i)*b1(i,1);
        si(1,i)=0;
    end
    for i=2:kluster
        aaaa(:,:)=a1(i-1,:,:);
        for j=1:st
            temp=0;
            for k=1:st
                temp(k)=del((i-1),k)*aaaa(k,j);
            end
            [temp1,temp2]=max(temp);
            del(i,j)=temp1*b1(j,i);
            si(i,j)=temp2;
        end
    end
    t=kluster;
    clear temp;
    temp=del(t,:);
    [temp1,temp2]=max(temp);
    psta=temp1;
    ssta(t)=temp2;
    for i=1:(kluster-1)
        t=t-1;
        ssta(t)=si((t+1),ssta(t+1));
    end
    sstar(l,:)=ssta(:);
end

%%%%%------------End of Viterbi-----------------%%%%%%%