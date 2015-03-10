phiqdata=xlsread('C:\Users\Baktha\Desktop\Personal\sampleData.xlsx');
k=1;
for i=1:2:143
    PHIp(:,k)=phiqdata(:,i);
    Qp(:,k)=phiqdata(:,i+1);
    Np(:,k)=phiqdata(:,i+1);
    k=k+1;
end

T=72;
S=50;
N=5;
tot=T*S; % total no of phi, q sets
m=0;


Qpmean=mean(Qp);

z=1;
for k=1:S
    for y=1:T
        %     if(mod(z,72)~=0)
        %         colmean=mod(y,72);
        %     else
        %         colmean=72;
        %     end
        colmean=y;
        if(y==1)
            if(Qp(k,y)<=Qpmean(colmean))
                state(z)=1;
            else
                state(z)=2;
            end
            z=z+1;
        elseif(y>1 && y<=9)
            if(Qp(k,y)<=Qpmean(colmean))
                state(z)=1;      %anomalies in PD region
            else
                state(z)=2;
            end
            z=z+1;
        elseif(y>9 && y<=18)    % upto 90 degrees
            if(Qp(k,y)<=Qpmean(colmean))
                state(z)=3;       %anomalies in PD region
            else
                state(z)=2;
            end
            z=z+1;
        elseif(y>18 && y<=27)    % beyond 90 deg
            if(Qp(k,y)<=Qpmean(colmean))
                state(z)=3;
            else
                state(z)=2;  % anomalies
            end
            z=z+1;
        elseif(y>27 && y<=36)     % upto 180 deg
            if(Qp(k,y)<=Qpmean(colmean))
                state(z)=3;
            else
                state(z)=4;   % anomalies
            end
            z=z+1;
        elseif(y>36 && y<=45) % beyond 180 deg
            if(Qp(k,y)<=Qpmean(colmean))
                state(z)=3;   % anomalies
            else
                state(z)=4;
            end
            z=z+1;
        elseif(y>45 && y<=54) % upto 270 deg
            if(Qp(k,y)<=Qpmean(colmean))
                state(z)=5;   % anomalies
            else
                state(z)=4;
            end
            z=z+1;
        elseif(y>54 && y<=63) % beyond 270
            if(Qp(k,y)<=Qpmean(colmean))
                state(z)=5;
            else
                state(z)=4; % anomalies
            end
            z=z+1;
        elseif(y>63 && y<=72) % upto 360 deg
            if(Qp(k,y)<=Qpmean(colmean))
                state(z)=5;
            else
                state(z)=4;
            end
            z=z+1;
        end
    end
end


k=1;
for x=1:S
    for y=1:T
        statexy(x,y)=state(k);
        k=k+1;
    end
end

% z=1;
% for k=1:S
% for y=1:tot
%     if(Qp(y)>0.1)
%         state(z)=1;
%         z=z+1;
%     elseif(y>=1 && y<=18)
%         state(z)=2;
%         z=z+1;
%     elseif(y>18 && y<=36)
%         state(z)=3;
%         z=z+1;
%     elseif(y>36 && y<=54)
%         state(z)=4;
%         z=z+1;
%     elseif(y>54 && y<=72)
%         state(z)=5;
%         z=z+1;
%     end
% end
% end


for i=1:N
    pii(i)=0;
    for m=1:tot
        if(state(m)==i)
            pii(i)=pii(i)+1;
        end
    end
    pi(i)=pii(i)/S;
end

for i=1:N
    for j=1:N
        aa(i,j)=0;
        for y=1:S
            for x=1:T-1
                if((state(y*x)==i)&&(statexy(y,(x+1))==j))
                    aa(i,j)=aa(i,j)+1;
                end
            end
        end
        a(i,j)=aa(i,j)/pii(i);
    end
end

%finding mean
k=0;
phitot=0;
qtot=0;
meanofstate=zeros(5,3);
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
   % meanofstate(i,:)=[phimean,qmean];
    ntot=5;
    nmean=ntot/k;
    meanofstate(i,:)= [phimean,qmean, nmean];
end

PQN=zeros(S,T,3);
m=1;
for i=1:S
    for j=1:T
        Qp(m)=Np(m);
        PQN(i,j,:)=[PHIp(m),Qp(m),Np(m)];
        m=m+1;
    end
end

%finding covariance
PQNtot=0;
t=0;
covar=zeros(5,3,3);
for m=1:N
    for i=1:S
        for j=1:T
            if(statexy(i,j)==m)
                PQNtot=PQNtot+PQN(m,:);
                tem=[PQN(i,j,1),PQN(i,j,2),PQN(i,j,3)];
                t=t+(((tem-meanofstate(m,:))')*(tem-meanofstate(m,:)));
                k=k+1;
            end
        end
    end
    covar(m,:,:)=t/k;
end

%finding b
M=3; % phi, q, n

for i=1:S
    tem=zeros(3,3);
    for j=1:N
        for m=1:T
            t1=1/((2*3.14)^(M/2));
            for x=1:3
                for y=1:3
                    tem(x,y)=covar(j,x,y);
                end
            end
            for x=1:3
                tem2(:,x)=PQN(i,m,x);
            end
            t2=1/((det(tem))^(0.5));
            t3=tem2-meanofstate(j,:);
            t4=inv(tem);
            t5=(tem2-meanofstate(j,:))';
            b(i,j,m)= t1*t2*exp(-0.5*t3*t4*t5);
        end
    end
end



