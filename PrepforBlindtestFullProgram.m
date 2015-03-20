clear all;
clc;
for iiii=1:9
    clc;
    clearvars -except iiii ;
    fileno=strcat('00',num2str(iiii));
    type='Surface';
    %%%%New1StartingCode;
    
    phiqdata=xlsread(strcat('J:\Datasets\',type,'\',type,'File',fileno,'-PQNTIdata.xlsx'));
    symb=21;
    Nsamp=size(phiqdata,1);
    counter=0;
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
    PHIpmod=PHIp;
    
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
    
    %%%%%%%UniformInitforLRmodel;
    
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
    
    
    
    clear b;
    b=B;
    a=A;
    N = st;
    K = symb;
    sum0=0;
    Ob=zeros(size(Observ(1,:)));
    T=length(Ob);
    Beta=zeros(T,N);
    Alpha=zeros(T,N);
    ZI=zeros(T,N,N);
    Gamma=zeros(T,N);
    E_T=zeros(1,N);
    E_I_J=zeros(N,N);
    E_Pi=zeros(1,N);
    E_A=zeros(N,N);
    N_E_A=zeros(N,N);
    E_B=zeros(N,K);
    N_E_B=zeros(N,K);
    sum1=zeros(K);
    p_v=zeros(N);
    status=zeros(1,N);
    nu=0.0;
    obiter=size(Ob,1);
    order=randperm(18,10);
    %order=randperm(18,18);
    order=[1:18];
    NoOfOb=1;
    looptimes=1;
    
    for xxx=1:NoOfOb
        Beta=zeros(T,N);
        Alpha=zeros(T,N);
        ZI=zeros(T,N,N);
        Gamma=zeros(T,N);
        E_T=zeros(1,N);
        E_I_J=zeros(N,N);
        E_Pi=zeros(1,N);
        E_A=zeros(N,N);
        E_B=zeros(N,K);
        Ob=Observ(order(xxx),:);
        SingleLoopCalcHMM;
        Alphastore(xxx,:,:)=Alpha;
        Betastore(xxx,:,:)=Beta;
        FinalProb(xxx)=ProbOgivenLfwd(end);
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
        N_E_B(N_E_B==0)=[0.00001];
        if sum(N_E_B==0)~=0
            for i=1:N
                coun(i)=sum(N_E_B(i,:)==0.00001);
            end
            for i=1:N
                N_E_B(i,:)=N_E_B(i,:)/(sum(E_B(i,:))+coun(i));
            end
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
    end
    mlogprob=zeros(1,looptimes);
    for j=1:looptimes
        for i=1:NoOfOb
            mlogprob(j)=mlogprob(j)+log(mlop(i,j))/ex;
        end
    end
    AijBjlDiffcalcMar16;
    
    filenamemodel=strcat('Lambda',type,'file',num2str(fileno),'.mat');
    save(filenamemodel,'Pi','N_E_A','N_E_B','New_A','New_B','mlop','mlogprob');
end