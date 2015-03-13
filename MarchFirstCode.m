clear all
clc

for mm=1:1
    clc;
    clearvars -except mm ;
    PI=pi;
    fileno=strcat('00',num2str(mm));
    type='Internal';
    %phiqdata=xlsread(strcat('J:\Universita di Bologna\Papers for modification in code\venkatesh sir 4 feb files\data of epoxy double v\epoxy void_ac test.',fileno,'PQNDATA.xlsx'));
    phiqdata=xlsread(strcat('J:\Datasets\',type,'\',type,'File',fileno,'-PQNdataCorrectedNQavg30equaldata.xlsx'));
    %phiqdata=xlsread(strcat('J:\Datasets\Interne\InternalFile',fileno,'-PQNdata.xlsx'));
    symb=21;
    Nsamp=30;
    counter=0;
    
    for j=1:Nsamp
        k=1;
        for i=1:3:214
            PHIp(j,k)=phiqdata(j,i);
            Qp(j,k)=phiqdata(j,i+1);
            if Qp(j,k)~=0
                counter=counter+1;
            end
            Np(j,k)=(phiqdata(j,i+2));
            k=k+1;
        end
    end
    Qpampl=abs(Qp);
    zz=k-1;
    qsum=0;
    
    % for j=1:Nsamp
    %     qsum=qsum+sum(Qpampl(j,:));
    % end
    
    
    wind=72;
    %qavg=qsum/counter;
    for j=1:Nsamp
        for k=1:zz
            if max(Qpampl(j,:))~=0
                Qpmod(j,k)=Qp(j,k)/max(Qpampl(:)); %changed this from Qpampl(j,:)
            else
                Qpmod(j,k)=0;
            end
        end
    end
    
    
    MinDistClusteringMar1
    ClusterandSeparatedataMar1
    CreateObservSeqMar1
    InitStateAssignMar1
    InitModelHMMmar1
%     StateSeqOptimMar1
%     HMMTrainingMar1
%     
%     save(strcat('var',type,'file',num2str(fileno),'.mat'), '');
%     delay(5);
end