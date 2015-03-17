clear all;
clc;
iiii=2;
fileno=strcat('00',num2str(iiii));
type='EpoxyInternal';
phiqdata=xlsread(strcat('J:\Datasets\',type,'\',type,'File',fileno,'-PQNTIdata.xlsx'));
%%%%%%%%%%%phiqdata=xlsread(strcat('J:\Datasets\',type,'\',type,'File',fileno,'-PQNdataCorrectedNQavg30equaldata.xlsx'));
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
        Np(j,k)=(phiqdata(j,i+2));
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


New2MinDistClusteringMar17
New3ClusterandSeparatedataMar17
New4CreateObservSeqMar17
New5InitStateAssignMar17
