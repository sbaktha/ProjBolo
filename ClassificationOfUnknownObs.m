fileno=strcat('00',num2str(iiii));
type='EpoxyInternal';
phiqdata=xlsread(strcat('J:\Datasets\',type,'\',type,'File',fileno,'-PQNTIdata.xlsx'));

New1StartingCode
New2MinDistClusteringMar17
New3ClusterandSeparatedataMar17
New4CreateObservSeqMar17

order=[1:18];
NoOfOb=18;
NoOfHmm=20;
for xxx=1:NoOfOb
    for yyy=1:NoOfHmm
        Beta=zeros(T,N);
        N = st;
        K = symb;
        T=length(Ob);
        Alpha=zeros(T,N);
        ZI=zeros(T,N,N);
        Gamma=zeros(T,N);
        E_T=zeros(1,N);
        E_I_J=zeros(1,N);
        E_Pi=zeros(1,N);
        E_A=zeros(N,N);
        E_B=zeros(N,K);
        filenamemodel=strcat('Lambda',type,'file',num2str(yyy),'.mat');
        load filenamemodel;
        Pi=E_Pi;
        a=New_A;
        b=New_B;
        SingleLoopCalcHMM;
        ProbModel(yyy)=ProbOgivenLfwd(end);
    end
    [ProbObs(xxx),Index(xxx)]=max(ProbModel(:));
end
count1=0;
count2=0;
count3=0;
for i=1:NoOfOb
    if Index(xxx)<=10
        DataType='Internal';
        count1=count1+1;
    elseif Index(xxx)>10 && Index(xxx)<=20
        DataType='External';
        count2=count2+1;
    elseif Index(xxx)>20 && Index(xxx)<=30
        DataType='Corona';
        count3=count3+1;
    end
end
