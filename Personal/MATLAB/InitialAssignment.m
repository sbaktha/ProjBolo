pqndata=pqndata1; %this returns phi, q, n for T different segments of S different samples
PHIp=pqndata(:,1);
Qp=pqndata(:,2);
Np=pqndata(:,3);

%N- No of states
%T- No of segments in one sample
%S- No of samples used for training

tot=T*S; % total no of phi, q, n sets
m=0; 
for y=1:S
    for x=1:T
        phi=PHIp(m);  % phi of the particular data set
        q=Qp(m);
        n=Np(m);    
        % define the states for each of the segments in all samples
        % should we use n or phi
        % state(y,x)=1-5;      
        %statetot(m)=1-5;
        m=m+1;
        % is it always a left-right model? or  can it be ergodic too? the state
        % changes from 5 to 1.. other than that it is left-right only.
    end
end

for i=1:N
    pii(i)=0;
    for m=1:tot
        if(statetot(m)==i)
         pii(i)=pii(i)+1;
        end
    end
    pi(i)=pii(i)/S;
end

for i=1:N
    for j=1:N
        aa(i,j)=0;
        for y=1:S
            for x=1:T
                if((state(y,x)==i)&&(state(y,(x+1))==j))
                    aa(i,j)=aa(i,j)+1;
                end
            end
        end
        a(i,j)=aa(i,j)/pii(i);
    end
end

%finding mean
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
    meanofstate(i)=[phimean,qmean,nmean];
end
%finding covariance
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

%finding b
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



        
                