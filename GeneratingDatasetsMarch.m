clear all;
clc;

for iii=3:3
    clearvars -except iii
    fileno=strcat('00',num2str(iii));
    type='Internal';
    filename1=strcat('J:\Datasets\',type,'\',type,'File',fileno);
    data=xlsread(strcat(filename1,'-PQNTIdata.xlsx'));
    siz=size(data);
    row=siz(1);
    col=siz(2);
    k=1;
    for i=1:4:72*4-3
        phi(:,k)=data(:,i);
        q(:,k)=data(:,i+1);
        no(:,k)=data(:,i+2);
        k=k+1;
    end
    
    m=1;
    if row<30
        for k=1:72
            qsu(m,k)=0;
            qav(m,k)=0;
            n(m,k)=0;
            nosu(m,k)=0;
            noav(m,k)=0;
            for j=i:size(data,1)
                qsu(m,k)=qsu(m,k)+q(j,k);
                nosu(m,k)=nosu(m,k)+no(j,k);
                if(q(j,k)~=0)
                    n(m,k)=n(m,k)+1;
                end
            end
            qav(m,k)=qsu(m,k)/size(data,1);
            noav(m,k)=nosu(m,k)/size(data,1);
        end
        m=m+1;
    else
        for i=1:floor(row/30):floor(row/30)*30
            for k=1:72
                qsu(m,k)=0;
                qav(m,k)=0;
                n(m,k)=0;
                nosu(m,k)=0;
                noav(m,k)=0;
                for j=i:i+floor(row/30)-1
                    qsu(m,k)=qsu(m,k)+(q(j,k));
                    nosu(m,k)=nosu(m,k)+no(j,k);
                    if(q(j,k)~=0)
                        n(m,k)=n(m,k)+1;
                    end
                end
                qav(m,k)=qsu(m,k)/ floor(row/30);
                noav(m,k)=nosu(m,k)/size(data,1);
            end
            m=m+1;
        end
        for i=floor(row/30)*30+1:row
            for k=1:72
                qsu(m,k)=0;
                qav(m,k)=0;
                n(m,k)=0;
                nosu(m,k)=0;
                noav(m,k)=0;
                qsu(m,k)=qsu(m,k)+(q(i,k));
                nosu(m,k)=nosu(m,k)+no(i,k);
                qav(m,k)=qsu(m,k)/ (row-floor(row/30)*30);
                noav(m,k)=nosu(m,k)/(row-floor(row/30)*30);
                m=m+1;
            end
        end
    end
    for i=1:size(qav,1) % for i=1:30
        for j=1:72
            combinedfinal(i,(3*j)-2)=phi(1,j);
            combinedfinal(i,(3*j)-1)=qav(i,j);
            combinedfinal(i,(3*j))=n(i,j);
        end
    end
    xlswrite(strcat(filename1,'-PQNdataCorrectedNQavg30equaldata.xlsx'),combinedfinal);
end