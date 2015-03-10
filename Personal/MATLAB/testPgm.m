T=36;
S=10;
N=5;
tot=T*S;
z=1;

for k=1:S
for y=1:tot
    if(y==1)
        state(z)=1;
        z=z+1;
    elseif(y>=1 && y<=10)
        state(z)=2;
        z=z+1;
    elseif(y>10 && y<=19)
        state(z)=3;
        z=z+1;
    elseif(y>18 && y<=29)
        state(z)=4;
        z=z+1;
    elseif(y>29 && y<=36)
        state(z)=5;
        z=z+1;
    end
end
end

for i=1:N
    pii(i)=0;
    for m=1:tot
        if(state(m)==i)
            pii(i)=pii(i)+1;
        end
    end
    pi(i)=pii(i)/S;
end

%
for i=1:N
    for j=1:N
        aa(i,j)=0;
        for m=1:tot-1
            if((state(m)==i)&&(state(m+1)==j))
                    aa(i,j)=aa(i,j)+1;
            end
        end
        a(i,j)=aa(i,j)/pii(i);
    end
end

%Finding Mean
k=0;
for i=1:N
    for m=1:tot
        if(state(m)==i)
            phitot=phitot+PHIp(m);
            qtot=qtot+Qp(m);
            ntot=ntot+Np(m);
            k=k+1;
        end
    end
    phimean=phitot/k;
    qmean=qtot/k;
    nmean=ntot/k;
    meanofstate(i)=[phimean;qmean;nmean];
end


for i=1:N
    for m=1:tot
        if(state(m)==i)
            phitot=phitot+PHIp(m);
            qtot=qtot+Qp(m);
            ntot=ntot+Np(m);
            k=k+1;
        end
    end
end



for j=1:N
    for m=1:tot
        t1=1/((2*pi)^(M/2));
        t2=1/((abs(covar(j)))^(0.5));
        t3=param(m)-meanofstate(j);
        t4=inv(covar(j));
        t5=(param(m)-meanofstate(j))^('t');
        b(j,m)= t1*t2*exp(-0.5*t3*t4*t5);
    end
end



