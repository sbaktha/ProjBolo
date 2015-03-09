clear all;
clc;
%aaaa=20;
%fileno=strcat('0',num2str(aaaa));
fileno='011';
%filename1=strcat('J:\Datasets\Interne\InternalFile',fileno);
filename1=strcat('J:\Datasets\data of epoxy double v\epoxy void_ac test.',fileno);
ampli=textread(strcat(filename1,'.pd2.A.txt'));
phase=textread(strcat(filename1,'.pd2.P.txt'));
times=textread(strcat(filename1,'.pd2.Ti.txt'));
tottime=floor(max(times));
for i=1:tottime+1
    ind=(times>(i-1) & times<=(i));
    last(i)=max(find(ind));
end
diff(1)=last(1);
for i=2:tottime+1
    diff(i)=(last(i)-last(i-1));
end
siz=max(diff);
firsttime=1;
inde=zeros(tottime+1,siz);
for i=1:tottime+1
    k=1;
    for j=firsttime:last(i)
        inde(i,k)=j;
        sec(i,k)=times(inde(i,k));
        secdata(i,k)=ampli(inde(i,k));
        k=k+1;
    end
    firsttime=last(i)+1;
end

qsorted=sort(abs(secdata),2);
m=1;
nonzeroqsorted={};
row=size(qsorted,1);
col=size(qsorted,2);
for i=1:row
    k=1;
    for j=1:col
        if qsorted(i,j)~=0
            temp(k)=qsorted(i,j);
            k=k+1;
        end
    end
    nonzeroqsorted{1,m}=[temp];
    m=m+1;
    clear temp
end

for i=1:row
    qnew=nonzeroqsorted{1,i};
    [c,d]=hist(qnew,unique(qnew));
    index10=find(c==10);
    if isempty(index10)~=1
        values10=d(index10);
        qnorm(i)=mean(values10);
    else
        qnorm(i)=1;  % Ask sir what should we do if there are no Q for which more than 10 pulses are there
    end
    clear qnew index10 values10 c d
end
k=1;
firsttime=1;
for i=1:row
    for j=firsttime:last(i)
        i;
        j;
        amplinorm(k,1)=ampli(k,1)/qnorm(i);
        k=k+1;
    end
    firsttime=last(i)+1;
end

