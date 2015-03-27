o(:,1)=[6.6096;-0.2507;0.3602];
o(:,2)=[2.3209;-1.6724;0.7802];
o(:,3)=[0.6991;-4.2347;0.9186];
o(:,4)=[0.2213;0.2285;0.9594];
o(:,5)=[0.2607;-2.2482;0.9866];
o(:,6)=[-0.0273;-2.3789;0.9985];
o(:,7)=[-0.0553;-1.8583;0.9938];
o(:,8)=[-0.1139;3.2653;0.9812];
o(:,9)=[0.1383;-1.0925;0.9860];
o(:,10)=[-0.0002;-0.0305;0.9919];
% replace the rs n clmns in kmeans
% write the indices in ascending order to revert
%O=[o1,o2,o3,o4,o5,o6,o7,o8,o9,o10];
flag=0;
cnt=0;
N=3;
T=6;
i=0;
Oblen=10;
c=zeros(3,N);
% count=zeros(1,10);
% for j=1:(Oblen-T)
%     for i=1:T
%         k(i)=sqrt((o(1,T+j)-o(1,i))^2+(o(2,T+j)-o(2,i))^2+(o(3,T+j)-o(3,i))^2;
%     end
%     [tem, index]=min(k);
%     count(T+j)=index;
% end
% for i=1:T
%     cnt=0;
%     sum=zeros(3,T);
%     sum(:,i)=0;
%     for j=T+1:Oblen
%         if(count(j)==i)
%             sum(:,i)=sum(:,i)+o(:,j);
%             cnt=cnt+1;
%         end
%     end
%     sum(:,i)=sum(:,i)+o(:,i);
%     c(:,i)=sum(:,i)/(cnt+1);
% end

for j=5:10
    for i=1:4
        k(i)=(sqrt((o(1,j)-o(1,i))^2+(o(2,j)-o(2,i))^2+(o(3,j)-o(3,i))^2));
    end
    [tem, index]=min(k);
    count(j)=index;
end
for j=5:10
    cnt=0;
    sum0=zeros(3,10);
    sum0(:,j)=0;
    for i=1:4
        if(count(j)==i)
            sum0(:,j)=sum0(:,j)+o(:,i);
            cnt=cnt+1;
        end
    end
    sum0(:,j)=sum0(:,j)+o(:,j);
    c(:,j)=sum0(:,j)/(cnt+1);
end
b=zeros(N,T);
for i=5:10
    b(:,i-4)=c(:,i);
end
clear c;
c=b;
for i=1:N
    st(:,i)=c(:,i);
end
count=0;
flag=0;
I=zeros(N,T);
cnt=0;
while(flag==0)
    stprev=st;
    d=zeros(N,T);
    for i=1:N
        dd=zeros(1,T);
        for j=1:T
            for k=1:3
                dd(1,j)=dd(1,j)+(st(k,i)-c(k,j))^2;
            end
            d(i,j)=sqrt(dd(1,j));
        end
    end
    
    Iprev=I;
    I=zeros(N,T);
    for i=1:T
        [temp,index]=min(d(:,i));
        I(index,i)=1;
    end
    if(Iprev==I)
        flag=1;
    end
    for i=1:N
        stt(:,i)=zeros(N,1);
        count=0;
        for j=1:T
            if(I(i,j)==1)
                stt(:,i)=stt(:,i)+c(:,j);
                count=count+1;
            end
        end
        st(:,i)=stt(:,i)/count;
    end
    %     if(stprev==st)
    %         flag=1;
    %     end
    cnt=cnt+1
end

for j=1:T
    for i=1:N
        if(I(i,j)==1)
            s1(j)=i;
        end
    end
end
    
