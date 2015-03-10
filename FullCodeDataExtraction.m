%-----TO READ DATA FROM EXCEL SHEET AND REARRANGE THEM INTO CYCLEWISE PHI
%AND Q VALUES------

fileno='001';
data=xlsread(strcat('J:\Datasets\Superficiali\SurfaceFile',fileno,'-rawPQdata.xlsx'));
%data=xlsread('J:\Universita di Bologna\sv sir files 28 jan\rawdataphiq8kvDIVepoxydoublevoidfile015.xlsx'); %Add appropriate path based on the file location
a=size(data);
m=a(1);
n=a(2); %Saves number of rows and columns as variables for future looping
app=round(data(1:m,1));%A matrix containing phi values rounded to nearest integer is created
q=zeros(1,360);
phi=[1:360];%Creates phi and q matrix for the first row
j=1;
for i=1:360
    if(j>1)&&(app(j,1)<app(j-1,1))
        break
    end
    if(phi(1,i)==app(j,1))
        phi(1,i)=data(j,1);
        q(1,i)=data(j,2);
        j=j+1;
    end
end
%Above code populates the first row for the phi and q matrices
d=j; %Place holder for value of j at the end of first cycle
phiappend=[1:360];
qappend=zeros(1,360);
%Additional row matrices created
%These row matrices will be successively appended to the first row at the
%end of each cycle
for count=j:m
    if(count>d)
        if(app(count,1)<app(count-1,1))
            phi=[phi;phiappend];
            q=[q;qappend];
            phiappend=[1:360];
            qappend=zeros(1,360);
        end
    end
    for i=1:360
        if(phiappend(1,i)==app(count,1))
            phiappend(1,i)=data(count,1);
            qappend(1,i)=data(count,2);
        end
    end
end
%All phi values are saved in 'phi' matrix and all q values are stored in
%'q' matrix to be accessed separately in order to ease the process of detecting
%peaks and differences
[x y]=size(phi);
y2=2*y;
combined=zeros(x,y2);
for i=1:x
    for j=1:y
        y2=2*y;
        combined(i,(2*j)-1)=phi(i,j);
        combined(i,(2*j))=q(i,j);
    end
end
%Above program combines the phi and q matrices to give one large matrix
%with phi and q values adjacent to each other.
%-----TO FIND MAX AND MINIMUM VALUES-----
n=input('Enter the desired size of window');
loop=n;
% gets the desired window size
B = zeros (1) % a matrix B which is 1x1
C = size(q,1) % C is a variable which is equal to number of rows in 'q'
phimax = zeros(C , 360/n); % matrix phimax with C rows and 360/n columns
phimin = zeros(C , 360/n);
qmax = zeros(C , 360/n);% matrix qmax with C rows and 360/n columns
qmin = zeros(C , 360/n);
qtempo = zeros (1,n); % a temporary matrix with 1 row and n columns
for k = 1:C % for loop for ensuring all the rows are covered
    p = 1; % variable to keep track of 'current column' in 'q'
    for j = 1:(360/n) % for loop for scanning all the columns in a row
        for i = 1:n
            qtempo (1,i) = q(k,p); % n elements taken in at a time from 'q'
            p = p+1;
        end
        [ maxi , locat] = max(qtempo(:)); % finding the largest element of the n taken elements
        B = maxi;
        qmax (k , j) = B; % storing that largest element in qmax matrix
        [row,column] = ind2sub(size(qtempo),locat);
        col = ( p - 1) - (n - column); % logic for finding the location of qmax(k,j)
        % in q so that element in phi with same coordinates
        % can be copied to phimax too
        phimax (k , j) = phi( k , col); % i. e for finding phimax corresponding to qmax.
        [ mini , locat1] = min(qtempo(:)); % for calculating qmin and phi min. same logic as above.
        B = mini;
        qmin (k , j) = B;
        [row1,column1] = ind2sub(size(qtempo),locat1);
        col1 = ( p - 1) - (n - column1);
        phimin ( k, j) = phi ( k , col1);
    end
end
%Below is the program for combining the max and min phi and q values obtained.
[x y]=size(phimax);
y2=2*y;
combinedmax=zeros(x,y2);%Zero matrices which will be populated with max and minimum phi and q values for the given window size
combinedmin=zeros(x,y2);
for i=1:x
    for j=1:y
        y2=2*y;
        combinedmax(i,(2*j)-1)=phimax(i,j);
        combinedmax(i,(2*j))=qmax(i,j);
        combinedmin(i,(2*j)-1)=phimin(i,j);
        combinedmin(i,(2*j))=qmin(i,j);
    end
end
%code for phi max etc. calc
deltat1mat= ones(C,360/n);
deltat1mat = (0.00001 * deltat1mat);
deltat2mat= ones(C,360/n);
deltat2mat = (0.00001 * deltat2mat);
numpulse = zeros(C , 360/n);
nat_d1=zeros(C,360/n);
nat_d2=zeros(C,360/n);
prep=zeros(C,360/n);
nrep=zeros(C,360/n);
remarks=zeros(C,360/n);

for k = 1:C % for loop for incrementing the rows in q matrix
    t1 = 0;
    t2 = 0;
    t3 = 0;
    q1 = 0;
    q2 = 0;
    q3 = 0;
    deltat1 = 0;
    deltat2 = 0;
    natural_deltat1 = 0; % extra column which displays +1 in case of homo. positive , -1 if vice versa , +0.5 in case of non-homo. positive to negative
    natural_deltat2 = 0; % and -0.5 vice versa and in case of no pulse or only 1 pulse exists, 0.
    for j = 1:360/n
        tempopos = zeros( 1 , n); % temporary matrices to hold value of i so that the corresponding column in phi matrix can be easily found
        % by using phi(k,tempopos(1,x))
        temponeg = zeros (1 , n);
        temposort = zeros (1, (2*n));% to sort the column values stored in tempopos and temponeg
        temposort2 = zeros (1, (2*n));
        posi = 0; % counters to keep count of number of positive and negative values.
        nega = 0; 
        x = 1; % represents the column of tempopos
        y= 1; % represents the column of temponeg
        for i = 1:n % for loop for incrementing columns in q matrix.
            B = q(k,( n*j) - (n -i)); % to check if the pulse is positive or negative
            if B > 0
                posi = posi + 1;
                tempopos(1,x) = (( n*j) - (n -i)); % the column value of the positive pulse is stored in tempopos
                x = x+1;
            elseif B < 0
                nega = nega + 1;
                temponeg(1,y) = (( n*j) - (n -i)); % the column value of the negative pulse is stored in temponeg
                y = y+1;
            else
                break;
            end
        end
        pospulse = posi % the total number of positive pulses in the window. this is for displaying.
        negpulse = nega
        numpuls = posi + nega 
        x = 1;
        y = 1;
        if ( posi ~= 0) && (nega == 0) % case:1. only positive pulses are present.
            if ( posi > 1) % checking if more than 1 positive pulse is available.
                t1 = phi(k , tempopos( 1 , x )); % t1 gets the value of phi corresponding to the 1st positive pulse. the column value of phi is obtained from tempopos
                q1 = q(k , tempopos( 1 , x ));
                x = x+1;
                posi = posi - 1;
                t2 = phi(k , tempopos( 1 , x ));% t2 gets the value of phi corresponding to the 2nd positive pulse. the column value of phi is obtained from tempopos
                q2 = q(k , tempopos( 1 , x ));
                x = x+1;
                posi = posi - 1;
                deltat1 = t2 - t1; % delta 1
                natural_deltat1 = +1;
                q3 = q2 - q1;
                if ( q3 < 0)
                    deltat1 = ((-1)*deltat1);
                end
                if posi == 0 % to check if 3rd pulse is available
                    deltat2 = 0;
                    break
                else
                    t3 = phi(k , tempopos( 1 , x ));% t3 gets the value of phi corresponding to the 3rd positive pulse. the column value of phi is obtained from tempopos
                    q3 = q(k , tempopos( 1 , x ));
                    posi = posi - 1;
                    deltat2 = t3 - t1; % delta 2
                    natural_deltat2 = +1;
                    q2 = q3 - q1;
                    if ( q2 < 0)
                        deltat2 = ((-1)*deltat2);
                    end
                end
            else % when only 1 positive pulse is available
            %    disp (' There are no homogeneous pulses!! ' );
                remarks(k,j)=1;
                break
            end
        elseif ( nega ~= 0) && (posi == 0) % case:2. only negative pulses are present
            if ( nega > 1) % same logic as prev if statement for positive pulse.
                t1 = phi(k , temponeg( 1 , y ));
                q1 = q(k , temponeg( 1 , y ));
                y = y+1;
                nega = nega - 1;
                t2 = phi(k , temponeg( 1 , y ));
                q2 = q(k , temponeg( 1 , y ));
                y = y+1;
                nega = nega - 1;
                deltat1 = t2 - t1;
                natural_deltat1 = -1;
                q3 = q2 - q1;
                if ( q3 < 0)
                    deltat1 = ((-1)*deltat1);
                end
                if nega == 0
                    deltat2 = 0;
                    continue
                else
                    t3 = phi(k , temponeg( 1 , y ));
                    q3 = q(k , temponeg( 1 , y ));
                    nega = nega - 1;
                    deltat2 = t3 - t1;
                    natural_deltat2 = -1;
                    q2 = q3 - q1;
                    if ( q2 < 0)
                        deltat2 = ((-1)*deltat2);
                    end
                end
            else
           %     disp (' There are no homogeneous pulses!! ' );
                remarks(k,j)=1;
                break;
            end
        elseif ( posi ~= 0) && (nega ~=0) % if both positive and negative pulses are available
            if ( posi == 1) && (nega == 1) % if only 1 positive and 1 negative pulse is available
                t1 = phi(k , tempopos( 1 , x ));
                t2 = phi(k , temponeg( 1 , y ));
                q1 = q(k , tempopos( 1 , x));
                q2 = q(k , temponeg( 1 , y ));
                deltat2 = 0;
                deltat1 = t2 - t1;
                natural_deltat1 = +0.5;
                if (deltat1 <0 )
                    deltat1 = (deltat1 * (-1));
                    natural_deltat1 = -0.5;
                end
                q3 = q2 - q1;
                if ( q3 < 0)
                    deltat1 = ((-1)*deltat1);
                end
          %      disp (' There are no homogeneous pulses!! ');
                remarks(k,j)=1;
                break
            else
                temposort = horzcat ( tempopos , temponeg);
                temposort2 = sort(temposort);
                temposort2(temposort2==0) = []; % to cut all the 0s which exist in the beg.
                %since the column values are arranged in ass. order thus leaving only the non-zero column values
                t1 = phi(k, temposort2(1,1));
                t2 = phi(k, temposort2(1,2));
                t3 = phi(k, temposort2(1,3));
                q1 = q(k, temposort2(1,1));
                q2 = q(k, temposort2(1,2));
                deltat1 = t2 - t1;
                q3 = q2 - q1;
                if ( q3 < 0)
                    deltat1 = ((-1)*deltat1);
                end
                deltat2 = t3 - t1;
                q3 = q(k, temposort2(1,3));
                q2 = q3 - q1;
                if ( q2 < 0)
                    deltat2 = ((-1)*deltat2);
                end
                t1 = q(k, temposort2(1,1));
                t2 = q(k, temposort2(1,2));
                t3 = q(k, temposort2(1,3));
                if ( t1 > 0)
                    if (t2 > 0)
                        natural_deltat1 = +1;
                    else
                        natural_deltat1 = +0.5;
                    end
                    if (t3 >0)
                        natural_deltat2 = +1;
                    else
                        natural_deltat2 = +0.5;
                    end
                else
                    if (t2 < 0)
                        natural_deltat1 = -1;
                    else
                        natural_deltat1 = -0.5;
                    end
                    if (t3 < 0)
                        natural_deltat2 = -1;
                    else
                        natural_deltat2 = -0.5;
                    end
                end
            end
        else
           % disp (' There are no pulses!!! ');
            remarks(k,j)=2;
        end
        deltat1mat(k,j)=deltat1;
        deltat2mat(k,j)=deltat2;
        prep(k,j)=pospulse;
        nrep(k,j)=negpulse;
        numpulse( k ,j) = numpuls;
        nat_d1(k,j)=natural_deltat1;
        nat_d2(k,j)=natural_deltat2;
    end
end
[x y]=size(phimax);
y2=11*y;
finalmatrix=zeros(x,y2);%Zero matrices which will be populated with max and minimum phi and q values for the given window size
for i=1:x
    for j=1:y
        y2=11*y;
        finalmatrix(i,(11*j)-10)=phimax(i,j);
        finalmatrix(i,(11*j)-9)=qmax(i,j);
        finalmatrix(i,(11*j)-8)=phimin(i,j);
        finalmatrix(i,(11*j)-7)=qmin(i,j);
        finalmatrix(i,(11*j)-6)=deltat1mat(i,j);
        finalmatrix(i,(11*j)-5)=deltat2mat(i,j);
        finalmatrix(i,(11*j)-4)=numpulse(i,j);
        finalmatrix(i,(11*j)-3)=nat_d1(i,j);
        finalmatrix(i,(11*j)-2)=nat_d2(i,j);
        finalmatrix(i,(11*j)-1)=prep(i,j);
        finalmatrix(i,(11*j))=nrep(i,j);
    end
end
disp(' Check finalmatrix for consolidated data on all windows' );
disp(' Check remarks matrix for remarks on each window.');
disp('2 means There were no pulses.');
disp('1 means There were no homogenous pulses.' );
