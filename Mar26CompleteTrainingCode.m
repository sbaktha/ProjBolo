clear all;
Observtot=[];
Statetotal=[];
for iiii=1:9
    clc;
    clearvars -except iiii Observtot Statetotal;
    if iiii<10
        fileno=strcat('00',num2str(iiii));
    else
        fileno=strcat('0',num2str(iiii));
    end
    type='EpoxyInternal';
    phiqdata=xlsread(strcat('J:\Datasets\',type,'\',type,'File',fileno,'-PQNTIdatacorrected.xlsx'));
    symb=21;
    Nsamp=size(phiqdata,1);
    counter=0;
    %New1StartingCode
    for j=1:Nsamp
        k=1;
        for i=1:4:288
            PHIp(j,k)=phiqdata(j,i);
            Qp(j,k)=phiqdata(j,i+1);
            if Qp(j,k)~=0
                counter=counter+1;
            end
            Np(j,k)=phiqdata(j,i+2);
            Ti(j,k)=phiqdata(j,i+3);
            k=k+1;
        end
    end
    Qpampl=abs(Qp);
    zz=k-1;
    qsum=0;
    wind=72;
    for j=1:Nsamp
        for k=1:zz
            if max(Qpampl(j,:))~=0
                Qpmod(j,k)=Qp(j,k)/max(Qpampl(:)); %changed this from Qpampl(j,:)
            else
                Qpmod(j,k)=0;
            end
        end
    end
    %%%%%%%New2MinDistClusteringMar17;
    z=0;
    ex = 18;
    PHIpmod=PHIp/360;
    
    for i=1:Nsamp
        kk=1;
        for j=1:72
            testdata(i,kk)=PHIpmod(i,j);
            testdata(i,kk+1)=Qpmod(i,j);
            kk=kk+2;
        end
    end
    
    samplesmat=sort(randperm(Nsamp,ex));
    for i=1:18
        w(i,:)=testdata(samplesmat(i),:);
    end
    
    w_o = zeros(size(w)); % create a zero matrix of same dimentions of weight
    dif = w - w_o; % initialise difference.
    count = 1;                   % initialise iteration count.
    while sum(dif(:)) ~= 0 && count ~= 500
        sum(dif(:));
        w_o = w; % remember the weights of previous iterations.
        for ii = 1:Nsamp
            for jj = 1:ex
                eq_dist(jj) = ((testdata(ii,:)-w(jj,:)) * ((testdata(ii,:)-w(jj,:))')) ; % equiledian distance
            end
            [temp,near_class(ii)] = min(eq_dist);  % find the cluster which is in minimum distance from the training exempler.
        end
        for ii = 1:ex
            [a,b] = find(near_class == ii);
            temp_sum = 0;
            for jj = b
                temp_sum = temp_sum + testdata(jj,:);
            end
            if sum(a) == 0
                count;
                ii;
            end
            w(ii,:) = temp_sum / sum(a);
        end
        dif = abs(w - w_o);
        count = count+1;
    end
    %%%%%%%%New3ClusterandSeparatedataMar17
    
    for i=1:ex
        for j=1:Nsamp
            eq_dist(j)=sqrt((w(i,:)-testdata(j,:))*((w(i,:)-testdata(j,:))'));
        end
        [temp,j1]=min(eq_dist);
        g(i,:)=testdata(j1,:);
    end
        
    counter1=0;
    for j=1:ex
        k=1;
        for i=1:2:143
            PHInew(j,k)=g(j,i);
            Qnew(j,k)=g(j,i+1);
            if Qnew(j,k)~=0
                counter1=counter1+1;
            end
            k=k+1;
        end
    end
    zz=k-1;
    
    for i=1:ex
        l=1;d=1;
        for j=1:144              % for j=1:216 if we use N also
            g1(i,l,d)=g(i,j);
            d=d+1;
            if d==3          % d==4 if we use N also
                d=1;
                l=l+1;
            end
        end
    end
    %%%%%%New4CreateObservSeqMar17
    Qamplnew=abs(Qnew);
    wind=72;
    for j=1:ex
        for k=1:zz
            if max(Qamplnew(j,:)~=0)
                Qpmodnew(j,k)=Qnew(j,k);
            else
                Qpmodnew(j,k)=0;
            end
        end
    end
    
    for j=1:ex
        for i=1:72
            zz=0;
            for k=0:0.05:1
                zz=zz+1;
                if abs(Qpmodnew(j,i))>=k
                    Observ(j,i)=zz;
                end
            end
        end
    end
    Observtot=[Observtot;Observ];
    %%%New5InitStateAssignMar1
    
    st=4;
    
    for k=1:ex
        clear g;
        g(:,:)=g1(k,:,:);
        for i=1:st
            c(i,:)=g(i,:);
        end
        cl=1;
        while cl==1
            cl=1;
            for i=1:st
                for j=1:wind
                    dm(i,j)=sqrt((c(i,:)-g(j,:))*((c(i,:)-g(j,:))'));
                end
            end
            sm=zeros(st,wind);
            [temp,temp1]=min(dm);
            for i=1:wind
                sm(temp1(i),i)=1;
            end
            tsm=sum(sm,2);
            for i=1:st
                if tsm(i,1)>1
                    temp2=0;
                    for j=1:wind
                        if sm(i,j)==1
                            temp2=temp2+g(j,:);
                        end
                    end
                    c(i,:)=temp2/tsm(i,1);
                end
            end
            
            for i=1:st
                for j=1:wind
                    dm(i,j)=sqrt((c(i,:)-g(j,:))*((c(i,:)-g(j,:))'));
                end
            end
            
            sm1=zeros(st,wind);
            [temp,temp1]=min(dm);
            for i=1:wind
                sm1(temp1(i),i)=1;
            end
            if sm==sm1
                cl=0;
            end
            
        end
        sumstate(:,k)=sum(sm1,2)
        ss(k,:)=temp1(1,:);
    end
    
    S=ss;
    Statetotal=[Statetotal;S];
end
details=strcat('Data',type,'ObservStatemod.mat');
save(details,'Observtot','Statetotal');

type='EpoxyInternal';
varname=strcat('Data',type,'ObservStatemod.mat');
load (varname);

%%%%%%%UniformInitforLRmodel;
st=4;
symb=21;
N=st;
K=symb;
Pi=zeros(1,N);
A=zeros(N);
B=zeros(N,K);

Pi(1)=1;
for i=2:N
    Pi(i)=0;
end

for i=1:N
    for j=1:N
        if j==i || j==i+1
            A(i,j)=0.5;
        end
    end
end

for i=1:N
    for j=1:K
        B(i,j)=1/K;
    end
end

%HmmCalc15
clear b;
b=B;
a=A;
N = st;
K = symb;
sum0=0;
Ob=zeros(size(Observtot(1,:)));
T=length(Ob);
sum1=zeros(K);
nu=0.0;
obiter=size(Ob,1);
%order=[1:18]
NoOfOb=size(Observtot,1);
order=randperm(NoOfOb,NoOfOb);
looptimes=25;
ex=NoOfOb;
for xxx=1:NoOfOb
    Beta=zeros(T,N);
    Alpha=zeros(T,N);
    ZI=zeros(T,N,N);
    Gamma=zeros(T,N);
    E_Pi=zeros(1,N);
    E_T=zeros(1,N);
    E_I_J=zeros(N,N);
    E_A=zeros(N,N);
    N_E_A=zeros(N,N);
    E_B=zeros(N,K);
    N_E_B=zeros(N,K);
    status=zeros(1,N);
    p_v=zeros(N);
    
    Ob=Observtot(order(xxx),:);
    
    %SingleLoopCalc
    looptimes=25;
    for yyy=1:looptimes
        %b
        mlop(xxx,yyy)=0;
        for i=1:N
            Alpha(1,i)=Pi(i) * b(i,Ob(1));
        end
        
        for t=1:T-1
            for j=1:N
                sum0=0;
                for i=1:N
                    sum0= sum0+ Alpha(t,i)*a(i,j);
                end
                Alpha(t+1,j)=sum0 * b(j,Ob(t+1));
            end
        end
        %Alpha;
        ProbOgivenLfwd(yyy)=0;
        for i=1:N
            ProbOgivenLfwd(yyy)=ProbOgivenLfwd(yyy)+Alpha(T,i);
        end
        
        for i=1:N
            Beta(T,i)=1;
        end
        
        for t=T-1:-1:1
            for i=1:N
                sum0=0;
                for j=1:N
                    sum0=sum0+(a(i,j)*Beta(t+1,j)*b(j,Ob(t+1)));
                end
                Beta(t,i)=sum0;
            end
        end
        %Beta;
        
        ProbOgivenLbkwd(yyy)=0;
        for i=1:N
            ProbOgivenLbkwd(yyy)=ProbOgivenLbkwd(yyy)+Pi(i)*b(i,Ob(1))*Beta(1,i);
        end
        
        
        ZI=zeros(T,N,N);
        kk=0;
        %Calculation of ZI values
        sum0=zeros(1,T);
        for t=1:T-1
            for i=1:N
                for j=1:N
                    sum0(t) = sum0(t) + (Alpha(t,i) *a(i,j) *b(j,Ob(t+1)) *Beta(t+1,j));
                end
            end
        end
        
        for t=1:T-1
            for i=1:N
                for j=1:N
                    nu=Alpha(t,i)*b(j,Ob(t+1))*Beta(t+1,j)*a(i,j);
                    ZI(t,i,j) = nu/sum0(t);
                end
            end
        end
        % Gamma alternate
        %     sum00=zeros(1,T);
        %     for t=1:T
        %         for i=1:N
        %             sum00(t)=sum00(t)+Alpha(t,i)*Beta(t,i);
        %         end
        %     end
        for t=1:T
            sum00=0;
            for i=1:N
                sum00=sum00+Alpha(t,i)*Beta(t,i);
            end
            for i=1:N
                Gammanew(t,i)=Alpha(t,i)*Beta(t,i)/sum00;
            end
        end
        
        
        for i=1:N
            Delta(1,i)=Pi(i)*b(i,Ob(1));
            Shi(1,i)=0;
        end
        for t=2:T
            for j=1:N
                for i=1:N
                    tofindmaxshi(i)=Delta(t-1,i)*a(i,j);
                    tofindmaxdel(i)=Delta(t-1,i)*a(i,j);
                end
                Delta(t,j)=max(tofindmaxdel)*b(j,Ob(t));
                [temp, Shi(t,j)]=max(tofindmaxshi);
                %                 Delta(t,j)=Delta(t,j)*b(j,Ob(t));
                clear tofindmax;
            end
        end
        for i=1:N
            tofindmax(i)=Delta(T,i);
        end
        [Pnew,StatenewViter(T)]=max(tofindmax);
        
        for t=T-1:-1:1
            StatenewViter(t)=Shi(t+1,StatenewViter(t+1));
        end
        
        for t=1:T-1
            for i=1:N
                sum0=0;
                for j=1:N
                    sum0 = sum0+ZI(t,i,j);
                end
                Gamma(t,i)=sum0;
            end
        end
        
        
        Gamma(T,:)=Gammanew(T,:);
        
        for i=1:N
            sum0=0;
            for t=1:T-1
                sum0= sum0 + Gamma(t,i);
            end
            E_T(i)=sum0;
        end
        
        for i=1:N
            for j=1:N
                sum0=0;
                for t=1:T-1
                    sum0= sum0+ZI(t,i,j);
                end
                E_I_J(i,j)=sum0;
            end
        end
        
        for i=1:N
            E_Pi(i)= Gamma(1,i);
        end
        for i=1:N
            for j=1:N
                sum0=0;
                nu=0;
                for t=1:T-1
                    sum0=sum0+ZI(t,i,j);
                    nu=nu+Gamma(t,i);
                end
                E_A(i,j) = (sum0 / nu)  ;
            end
        end
        for j=1:N     % number of states
            sum2=0;
            for kk=1:K
                sum1(kk)=0;
            end
            for t=1:T  %to traverse the ob  servation sequence...
                sum2 = sum2+ Gamma(t,j);
                for kk=1:K
                    if(Ob(t) == kk)
                        sum1(kk)= sum1(kk) + Gamma(t,j);
                        break;
                    end
                end
            end
            for kk=1:K
                E_B(j,kk) = (sum1(kk))/sum2;
            end
        end
        
        mlop(xxx,yyy)=mlop(xxx,yyy)+ProbOgivenLbkwd(yyy);
        
        a=E_A;
        b=E_B;
        b(E_B<0.0001)=0.0001;
        for i=1:4
            b(i,:)=bsxfun(@rdivide,b(i,:),sum(b(i,:)));
        end
    end
    %mlop(xxx,:);
    
    Alphastore(xxx,:,:)=Alpha;
    Betastore(xxx,:,:)=Beta;
    FinalProb(xxx)=ProbOgivenLfwd(end);
    %mlogprob(xxx)=sum(log(mlop(xxx,:))/obiter);
    %Normalization
    sum2 = 0;
    sum3 = 0;
    pp = 0;
    pp1 = 0;
    for i=1:N
        if(i==status(pp+1)) % status(pp)
            pp=pp+1;
            for j=1:N
                N_E_A(i,j)=E_A(i,j);
            end
        else
            sum3=0;
            sum2=0;
            pp1=0;
            for j=1:N
                if(j==status(pp1+1)) %status(pp1)
                    pp1=pp1+1;
                    sum3= sum3 + a(i,j);
                else
                    sum2= sum2 + E_A(i,j);
                end
            end
            pp1=0;
            for j=1:N
                if(j~=status(pp1+1))  %status(pp1)
                    N_E_A(i,j)=(1-sum3)*(E_A(i,j)/sum2);
                else
                    pp1=pp1+1;
                    N_E_A(i,j)=a(i,j);
                end
            end
        end
    end
    
    
    N_E_B=E_B;
    N_E_B(N_E_B<0.0001)=0.0001;
    for i=1:4
        N_E_B(i,:)=bsxfun(@rdivide,N_E_B(i,:),sum(N_E_B(i,:)));
    end
    a=N_E_A;
    b=N_E_B;
    Pi=E_Pi;
    
    disp('After Normalization:');
    fprintf('\n');
    disp('The estimated state transition matrix is:');
    for i=1:N
        for j=1:N
            fprintf('%.8f',N_E_A(i,j));
            fprintf('    ');
        end
        fprintf('\n');
    end
    fprintf('\n');
    disp('The estimated probability matrix is:');
    for i=1:N
        for j=1:K
            fprintf('%.8f',N_E_B(i,j));
            fprintf('    ');
        end
        fprintf('\n');
    end
    astore(xxx,:,:)=E_A;
    bstore(xxx,:,:)=E_B;
    anstore(xxx,:,:)=N_E_A;
    bnstore(xxx,:,:)=N_E_B;
end

mlogprob=zeros(1,looptimes);
for j=1:looptimes
    for i=1:NoOfOb
        mlogprob(j)=mlogprob(j)+log(mlop(i,j))/ex;
    end
end
AijBjlDiffcalcMar23;



