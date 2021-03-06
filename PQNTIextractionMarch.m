%-----TO READ DATA FROM EXCEL SHEET AND REARRANGE THEM INTO CYCLEWISE PHI
%AND Q VALUES------
clear all;
clc;
for aaaa=10:20
    clearvars -except aaaa
    fileno=strcat('0',num2str(aaaa));
    type='Internal';
    filename1=strcat('J:\Datasets\',type,'\',type,'File',fileno);
    data=xlsread(strcat(filename1,'-rawPQTIdata.xlsx'));
    siz=size(data);
    m=siz(1);
    n=siz(2);               %Saves number of rows and columns as variables for future looping
    app=round(data(1:m,1)); %A matrix containing phi values rounded to nearest integer is created
    phi=[1:360];            %Creates phi, q, n, times matrix for the first row
    q=zeros(1,360);
    pti=zeros(1,360);
    no=zeros(1,360);
    j=1;
    clear phi; %
    phi=[1:360]; %
    flag=0;
    i=1;
    while i<=(360+flag)
        x=i-flag;
        if(j>1)&&(data(j,1)<data(j-1,1))
            break;
        end
        if(phi(1,x)==app(j,1))
            phi(1,x)=data(j,1);
            q(1,x)=data(j,2);
            pti(1,x)=data(j,3);
            j=j+1;
            if(round(phi(1,x))==app(j,1))
                phi(1,x-1)=data(j-1,1);
                q(1,x-1)=data(j-1,2);
                pti(1,x-1)=data(j-1,3);
                phi(1,x)=round(phi(1,x));
                flag=flag+1;
            end
        end
        i=i+1;
    end
    %Above code populates the first row for the phi and q matrices
    
    d=j; %Place holder for value of j at the end of first cycle
    phiappend=[1:360];
    qappend=zeros(1,360);
    ptiappend=zeros(1,360);
    
    %Additional row matrices created
    %These row matrices will be successively appended to the first row at the
    %end of each cycle
    
    for count=j:m
        if(count>d)
            if(app(count,1)<app(count-1,1))
                phi=[phi;phiappend];
                q=[q;qappend];
                pti=[pti;ptiappend];
                phiappend=[1:360];
                qappend=zeros(1,360);
                ptiappend=zeros(1,360);
            end
        end
        for i=1:360
            if(phiappend(1,i)==app(count,1))
                phiappend(1,i)=data(count,1);
                qappend(1,i)=data(count,2);
                ptiappend(1,i)=data(count,3);
            end
        end
    end
    npulf=(q~=0);
    
    %All phi values are saved in 'phi' matrix and all q values are stored in
    %'q' matrix to be accessed separately in order to ease the process of detecting
    %peaks and differences
    [x y]=size(phi);
    y2=4*y;
    combined=zeros(x,y2);
    for i=1:x
        for j=1:y
            y2=4*y;
            combined(i,(4*j)-3)=phi(i,j);
            combined(i,(4*j)-2)=q(i,j);
            combined(i,(4*j)-1)=npulf(i,j);
            combined(i,(4*j))=pti(i,j);
        end
    end
    % Above program combines the phi and q matrices to give one large matrix
    % with phi and q values adjacent to each other.
    
    %-----TO FIND MAX AND MINIMUM VALUES-----
    
    n=5;
    loop=n;
    B = zeros (1);              % a matrix B which is 1x1
    C = size(q,1);              % C is a variable which is equal to number of rows in 'q'
    phimax = zeros(C , 360/n);  % matrix phimax with C rows and 360/n columns
    phimin = zeros(C , 360/n);
    qmax = zeros(C , 360/n);    % matrix qmax with C rows and 360/n columns
    qmin = zeros(C , 360/n);
    npulses=zeros(C,360/n);
    qtempo = zeros (1,n);       % a temporary matrix with 1 row and n columns
    for k = 1:C                 % for loop for ensuring all the rows are covered
        p = 1;                  % variable to keep track of 'current column' in 'q'
        for j = 1:(360/n)       % for loop for scanning all the columns in a row
            countn=0;
            for i = 1:n
                qtempo (1,i) = q(k,p);          % n elements taken in at a time from 'q'
                countn=countn+npulf(k,p);
                p = p+1;
            end
            [ maxi , locat] = max(qtempo(:));   % finding the largest element of the n taken elements
            B = maxi;
            qmax (k , j) = B;                   % storing that largest element in qmax matrix
            [row,column] = ind2sub(size(qtempo),locat);
            col = ( p - 1) - (n - column);      % logic for finding the location of qmax(k,j)
            % in q so that element in phi with same coordinates
            % can be copied to phimax too
            phimax (k , j) = phi( k , col);     % i. e for finding phimax corresponding to qmax.
            ptimes(k,j)=pti(k,col);
            npulses(k,j)=countn;
            [ mini , locat1] = min(qtempo(:));  % for calculating qmin and phi min. same logic as above.
            B = mini;
            qmin (k , j) = B;
            ptimes1(k,j)=pti(k,locat1);
            [row1,column1] = ind2sub(size(qtempo),locat1);
            col1 = ( p - 1) - (n - column1);
            phimin ( k, j) = phi ( k , col1);
            ptimes1(k,j)=pti(k,col1);
            if abs(qmin(k,j))> abs(qmax(k,j))
                qfinal(k,j)=qmin(k,j);
                ptifinal(k,j)=ptimes1(k,j);
            else
                qfinal(k,j)=qmax(k,j);
                ptifinal(k,j)=ptimes(k,j);
            end
        end
    end
    
    %Below is the program for combining the max and min phi and q values obtained.
    [x y]=size(phimax);
    combinedfinal=zeros(x,4*x);
    for i=1:x
        for j=1:y
            combinedfinal(i,(4*j)-3)=phimax(i,j);
            combinedfinal(i,(4*j)-2)=qfinal(i,j);
            combinedfinal(i,(4*j)-1)=npulses(i,j);
            combinedfinal(i,(4*j))=ptifinal(i,j);
        end
    end
    
    [Nsamp,wind]=size(combinedfinal);
    m=0;
    a=1;b=1;
    for j=1:Nsamp
        k=1;
        for i=1:4:(wind-3)
            Qp(j,k)=combinedfinal(j,i+1);
            if(Qp(j,k)~=0)
                qnew(1,a)=Qp(j,k);
                a=a+1;
                % b=b+1;
            end
            k=k+1;
        end
        m=m+k;
    end
    
    [c,d]=hist(qnew,unique(qnew));
    nop=[d;c];
    
    y=1;
    for i=1:72
        PQNTI(:,y)=phimax(:,i);
        y=y+1;
        PQNTI(:,y)=qfinal(:,i);
        y=y+1;
        PQNTI(:,y)=npulses(:,i);
        y=y+1;
        PQNTI(:,y)=ptifinal(:,i);
        y=y+1;
    end
    
    xlswrite(strcat(filename1,'-PQNTIdata.xlsx'),PQNTI);
end