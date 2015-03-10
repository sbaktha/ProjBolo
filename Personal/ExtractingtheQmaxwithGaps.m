clear all;
clc;
aaaa=4;
fileno=strcat('00',num2str(aaaa));
filename1=strcat('J:\Datasets\Interne\InternalFile',fileno);
ampli=textread(strcat(filename1,'.pdb.A.txt'));
phase=textread(strcat(filename1,'.pdb.P.txt'));
times=textread(strcat(filename1,'.pdb.Ti.txt'));
tottime=floor(max(times));
loop=tottime+1;

for i=1:loop
    ind=(times>=(i-1) & times<=(i));
    if sum(ind)~=0
        last(i)=max(find(ind));
        first(i)=min(find(ind));
    end
end
diff(1)=last(1);
for i=2:loop %tottime+1 inst of loop
    diff(i)=(last(i)-last(i-1));
end
siz=max(diff);
firsttime=1;
inde=zeros(loop,siz);
for i=1:loop
    k=1;
    for j=first(i):last(i)
        inde(i,k)=j;
        if inde(i,k)~=0
            sec(i,k)=times(inde(i,k));
            secdata(i,k)=ampli(inde(i,k));
        else
            sec(i,k)=0;
            secdata(i,k)=0;
        end
            k=k+1;
    end
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
            end
            k=k+1;
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
        for j=first(i):last(i)
            if j~=0
            amplinorm(k,1)=ampli(k,1)/qnorm(i);
            k=k+1;
            end
        end
    end
    