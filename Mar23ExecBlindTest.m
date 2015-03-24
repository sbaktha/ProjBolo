clear all;
clc;

iiii=1;
fileno=strcat('00',num2str(iiii));
type='Corona';
phiqdata=xlsread(strcat('J:\Datasets\',type,'\',type,'File',fileno,'-PQNTIdata.xlsx'));
Class=['Internal';'Surfacee';'Coronaaa'];
New1StartingCode
New2MinDistClusteringMar17
New3ClusterandSeparatedataMar17
New4CreateObservSeqMar17

order=[1:18];
NoOfOb=18;
NoOfHmm=3;
N = st;
K = symb;
T=size(Observ(1,:),2);
for xxx=1:NoOfOb
    Ob=Observ(order(xxx),:);
    for yyy=1:NoOfHmm
        Beta=zeros(T,N);
        Alpha=zeros(T,N);
        ZI=zeros(T,N,N);
        Gamma=zeros(T,N);
        E_T=zeros(1,N);
        E_I_J=zeros(N,N);
        E_Pi=zeros(1,N);
        E_A=zeros(N,N);
        E_B=zeros(N,K);
        if yyy==1
            type='Internal';
        elseif yyy==2
            type='Surface';
        elseif yyy==3
            type='Corona';
        end
        filenamemodel=strcat('Lambda',type,'file',num2str(fileno),'.mat');
        load(filenamemodel);
        Pi=[1,0,0,0];
        a=New_A;
        b=New_B;
        b(b==0)=0.0001;
        for i=1:4
            b(i,:)=bsxfun(@rdivide,b(i,:),sum(b(i,:)));
        end
        
        %SingleLoopCalcHMM;
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
        
        ProbModel(yyy)=0;
        for i=1:N
            ProbModel(yyy)=ProbModel(yyy)+Alpha(T,i);
        end
    end
    [ProbObs(xxx),Index(xxx)]=max(ProbModel(:));
end
count1=0;
count2=0;
count3=0;
for xxx=1:NoOfOb
    if Index(xxx)==1
        count1=count1+1;
    elseif Index(xxx)==2
        count2=count2+1;
    elseif Index(xxx)==3
        count3=count3+1;
    end
end

Ans=[count1,count2,count3];
[prob,index]=max(Ans)
Class(index,:)