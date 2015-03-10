clear all;
clc;
fileno='015';
%data=xlsread(strcat('J:\Universita di Bologna\Papers for modification in code\venkatesh sir 4 feb files\data of epoxy double v\epoxy void_ac test.',(fileno),'PQNDATA.xlsx'));
data=xlsread(strcat('J:\Datasets\data of epoxy double v\epoxy void_ac test.',fileno,'PQNDATA.xlsx'));
[Nsamp,wind]=size(data);
m=0;
a=1;b=1;
for j=1:Nsamp
    k=1;
    for i=1:3:(wind-2)
        phi(j,k)=data(j,i);
        q(j,k)=data(j,i+1);
      %  Np(j,k)=data(j,i+2);
        if(q(j,k)~=0)
            qnew(1,a)=q(j,k);
            a=a+1;
           % b=b+1;
        end
        k=k+1;
    end
    m=m+k;
end

[c,d]=hist(qnew,unique(qnew));
nop=[d;c];  
for i=1:size(q,1)
    for j=1:size(q,2)
        for k=1:size(nop,2)
            if q(i,j)==nop(1,k)
                Np(i,j)=nop(2,k);
            end
        end
    end
end

 t=0:1:360;
 plot(t,5*sind(t))
 hold on
 for i=1:3373
     scatter(phi(i,:),q(i,:),3,'k.')
 end