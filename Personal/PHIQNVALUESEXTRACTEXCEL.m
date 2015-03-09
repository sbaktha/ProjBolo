%-----TO READ DATA FROM EXCEL SHEET AND REARRANGE THEM INTO CYCLEWISE PHI
%AND Q VALUES------
clear all;
for aaaa=1:9
    clearvars -except aaaa
    fileno=strcat('00',num2str(aaaa));
    ampli=textread(strcat('J:\Datasets\Corona\CoronaFile',fileno,'.pdb.A.txt'));
    phase=textread(strcat('J:\Datasets\Corona\CoronaFile',fileno,'.pdb.P.txt'));
    times=textread(strcat('J:\Datasets\Corona\CoronaFile',fileno,'.pdb.Ti.txt'));
    %times=textread('J:\Universita di Bologna\Papers for modification in code\venkatesh sir 4 feb files\data of epoxy double v\epoxy void_ac test.024.pd2.Ti.txt');
    var(:,2)=ampli;
    var(:,1)=phase;
    var(:,3)=phase;
    xlswrite(strcat('J:\Datasets\Corona\CoronaFile',fileno,'-rawPQdata.xlsx'),var);
    
    %clear all;
    %clc;
    %fileno='015';
    %data=xlsread(strcat('J:\Universita di Bologna\Papers for modification in code\venkatesh sir 4 feb files\data of epoxy double v\epoxy void_ac test.',(fileno),'.xlsx')); %Add appropriate path based on the file location
    %data=xlsread('J:\Universita di Bologna\sv sir files 28 jan\rawdataphiq8kvDIVepoxydoublevoidfile015.xlsx');
    filename1=strcat('J:\Datasets\Corona\CoronaFile',fileno);
    data=xlsread(strcat(filename1,'-rawPQdata.xlsx'));
    siz=size(data);
    m=siz(1);
    n=siz(2);%Saves number of rows and columns as variables for future looping
    app=round(data(1:m,1));%A matrix containing phi values rounded to nearest integer is created
    %no=zeros(1,360);
    phi=[1:360];%Creates phi and q matrix for the first row
    j=1;
    q=zeros(1,360);
    for i=1:360
        if(j>1)&&(app(j,1)<app(j-1,1))
            break;
        end
        if(phi(1,i)==app(j,1))
            phi(1,i)=data(j,1);
            q(1,i)=data(j,2);
            %no(1,i)=data(j,3);
            j=j+1;
        end
    end
    %Above code populates the first row for the phi and q matrices
    d=j; %Place holder for value of j at the end of first cycle
    phiappend=[1:360];
    qappend=zeros(1,360);
    %noappend=zeros(1,360);
    %Additional row matrices created
    %These row matrices will be successively appended to the first row at the
    %end of each cycle
    for count=j:m
        if(count>d)
            if(app(count,1)<app(count-1,1))
                phi=[phi;phiappend];
                q=[q;qappend];
                %  no=[no;noappend];
                phiappend=[1:360];
                qappend=zeros(1,360);
                %   noappend=zeros(1,360);
            end
        end
        for i=1:360
            if(phiappend(1,i)==app(count,1))
                phiappend(1,i)=data(count,1);
                qappend(1,i)=data(count,2);
                %           noappend(1,i)=data(count,3);
            end
        end
    end
    
    %All phi values are saved in 'phi' matrix and all q values are stored in
    %'q' matrix to be accessed separately in order to ease the process of detecting
    %peaks and differences
    [x y]=size(phi);
    y2=3*y;
    combined=zeros(x,y2);
    for i=1:x
        for j=1:y
            y2=3*y;
            combined(i,(3*j)-2)=phi(i,j);
            combined(i,(3*j)-1)=q(i,j);
            %combined(i,(3*j))=Np(i,j);
        end
    end
    %Above program combines the phi and q matrices to give one large matrix
    %with phi and q values adjacent to each other.
    
    
    %-----TO FIND MAX AND MINIMUM VALUES-----
    %n=input('Enter the desired size of window');
    n=5;
    loop=n;
    % gets the desired window size
    B = zeros (1); % a matrix B which is 1x1
    C = size(q,1); % C is a variable which is equal to number of rows in 'q'
    phimax = zeros(C , 360/n); % matrix phimax with C rows and 360/n columns
    phimin = zeros(C , 360/n);
    qmax = zeros(C , 360/n);% matrix qmax with C rows and 360/n columns
    qmin = zeros(C , 360/n);
    %nomax=zeros(C,360/n);
    %nomin=zeros(C,360/n);
    qtempo = zeros (1,n); % a temporary matrix with 1 row and n columns
    for k = 1:C         % for loop for ensuring all the rows are covered
        p = 1;              % variable to keep track of 'current column' in 'q'
        for j = 1:(360/n) % for loop for scanning all the columns in a row
            for i = 1:n
                qtempo (1,i) = q(k,p); % n elements taken in at a time from 'q'
                p = p+1;
            end
            [ maxi , locat] = max(qtempo(:)); % finding the largest element of the n taken elements
            B = maxi;
            qmax (k , j) = B;     % storing that largest element in qmax matrix
            [row,column] = ind2sub(size(qtempo),locat);
            col = ( p - 1) - (n - column); % logic for finding the location of qmax(k,j)
            % in q so that element in phi with same coordinates
            % can be copied to phimax too
            phimax (k , j) = phi( k , col); % i. e for finding phimax corresponding to qmax.
            %       nomax(k,j)= no(k,col);
            [ mini , locat1] = min(qtempo(:)); % for calculating qmin and phi min. same logic as above.
            B = mini;
            qmin (k , j) = B;
            [row1,column1] = ind2sub(size(qtempo),locat1);
            col1 = ( p - 1) - (n - column1);
            phimin ( k, j) = phi ( k , col1);
            %      nomin(k,j)=no(k,col1);
            if abs(qmin(k,j))> abs(qmax(k,j))
                qfinal(k,j)=qmin(k,j);
            else
                qfinal(k,j)=qmax(k,j);
            end
        end
    end
    %Below is the program for combining the max and min phi and q values obtained.
    [x y]=size(phimax);
    y2=3*y;
    %combinedmax=zeros(x,y2);%Zero matrices which will be populated with max and minimum phi and q values for the given window size
    %combinedmin=zeros(x,y2);
    combinedfinal=zeros(x,y2);
    for i=1:x
        for j=1:y
            y2=3*y;
            % combinedmax(i,(3*j)-2)=phimax(i,j);
            % combinedmax(i,(3*j)-1)=qmax(i,j);
            %combinedmax(i,(3*j))=Np(i,j);
            % combinedmin(i,(3*j)-2)=phimin(i,j);
            % combinedmin(i,(3*j)-1)=qmin(i,j);
            %combinedmin(i,(3*j))=Np(i,j);
            combinedfinal(i,(3*j)-2)=phimax(i,j);
            combinedfinal(i,(3*j)-1)=qfinal(i,j);
            %combinedfinal(i,(3*j))=Np(i,j);
        end
    end
    
    %xlswrite(strcat('J:\Universita di Bologna\Papers for modification in code\venkatesh sir 4 feb files\data of epoxy double v\epoxy void_ac test.',(fileno),'PQNDATA.xlsx'),combinedfinal);
    
    %data=xlsread(strcat('J:\Universita di Bologna\Papers for modification in code\venkatesh sir 4 feb files\data of epoxy double v\epoxy void_ac test.',(fi leno),'PQNDATA.xlsx'));
    [Nsamp,wind]=size(combinedfinal);
    m=0;
    a=1;b=1;
    for j=1:Nsamp
        k=1;
        for i=1:3:(wind-2)
            % phi(j,k)=data(j,i);
            Qp(j,k)=combinedfinal(j,i+1);
            %  Np(j,k)=data(j,i+2);
            if(Qp(j,k)~=0)
                qnew(1,a)=Qp(j,k);
                a=a+1;
                % b=b+1;
            end
            k=k+1;
        end
        m=m+k;
    end
    Np=zeros(size(Qp));
    [c,d]=hist(qnew,unique(qnew));
    nop=[d;c];
    for i=1:size(Qp,1)
        for j=1:size(Qp,2)
            for k=1:size(nop,2)
                if Qp(i,j)==nop(1,k)
                    Np(i,j)=nop(2,k);
                end
            end
        end
    end
    
    y=1;
    for i=1:72
        PQN(:,y)=phimax(:,i);
        y=y+1;
        PQN(:,y)=Qp(:,i);
        y=y+1;
        PQN(:,y)=Np(:,i);
        y=y+1;
    end
    %xlswrite(strcat('J:\Universita di Bologna\Papers for modification in code\venkatesh sir 4 feb files\data of epoxy double v\epoxy void_ac test.',(fileno),'PQNDATA.xlsx'),PQN);
    xlswrite(strcat(filename1,'-PQNdata.xlsx'),PQN);
end