clear all;
clc;
for aaaa=1:1
    clearvars -except aaaa
    aaaa=1;
    fileno=strcat('00',num2str(aaaa));
    type='Internal';
    filename1=strcat('J:\Datasets\',type,'\',type,'File',fileno);
    ampli=textread(strcat(filename1,'.pdb.A.txt'));
    phase=textread(strcat(filename1,'.pdb.P.txt'));
    times=textread(strcat(filename1,'.pdb.Tioffset.txt'));
    var(:,1)=phase;
    var(:,2)=ampli;
    var(:,3)=times;
    xlswrite(strcat(filename1,'-rawPQdata.xlsx'),var);
    data=xlsread(strcat(filename1,'-rawPQdata.xlsx'));
    siz=size(data);
    m=siz(1);
    n=siz(2);                   %Saves number of rows and columns as variables for future looping
    app=round(data(1:m,1));     %A matrix containing phi values rounded to nearest integer is create
    j=1;
    
    % % Using phase
    mm=1;
    noOfpti(mm)=0;
    k=1;
    while  k<size(phase,1)
        if phase(k+1)>phase(k)
            noOfpti(mm)=noOfpti(mm)+1;
            k=k+1;
        else
            rowlimit=(0.02*(1+phase(k)/360));
            if(times(k+1)-times(k))<=rowlimit
                mm=mm+1;
                noOfpti(mm)=0;
            elseif(times(k+1)-times(k))>rowlimit
                rowlimit=rowlimit-(0.02*(360-phase(k))/360);
                while(rowlimit>0)
                    mm=mm+1;
                    noOfpti(mm)=0;
                    rowlimit=rowlimit-0.02;
                end
            end
            k=k+1;
        end
    end
        noofrows=mm;
        
        phi=[1:360];
        phiappend=[1:360];
        q=zeros(noofrows,360);
        pti=zeros(noofrows,360);
        
        for k=1:noofrows-1
            phi=[phi;phiappend];
        end
        
            mm=1;
            totpt=0;
        for k=1:noofrows
            totpt=totpt+noOfpti(k);
            for i=1:360
                if(mm<=totpt)
                    if(phi(k,i)==app(mm,1))
                        phi(k,i)=data(mm,1);
                        q(k,i)=data(mm,2);
                        pti(k,i)=data(mm,3);
                        mm=mm+1
                    end
                end
            end
        end
        npulf=(q~=0);
        
        
        % All phi values are saved in 'phi' matrix and all q values are stored in
        % 'q' matrix to be accessed separately in order to ease the process of detecting
        % peaks and differences
        [x y]=size(phi);
        y2=3*y;
        combined=zeros(x,y2);
        for i=1:x
            for j=1:y
                y2=3*y;
                combined(i,(4*j)-3)=phi(i,j);
                combined(i,(4*j)-2)=q(i,j);
                combined(i,(4*j)-1)=npulf(i,j);
                combined(i,(4*j))=pti(i,j);
            end
        end
        %Above program combines the phi and q matrices to give one large matrix
        %with phi and q values adjacent to each other.
        
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
        
        xlswrite(strcat(filename1,'-PQNTIdatacorrect.xlsx'),PQNTI);
    end