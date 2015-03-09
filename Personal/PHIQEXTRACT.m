
% code for phi max etc. calc
deltat1mat= ones(C,360/n);
deltat1mat = (0.00001 * deltat1mat)
deltat2mat= ones(C,360/n);
deltat2mat = (0.00001 * deltat2mat)
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
    natural_deltat1 = 0; % extra column which displays +! in case of homo. positive , -1 if vice versa , +0.5 in case of non-homo. positive to negative
    natural_deltat2 = 0; %and -0.5 vice versa and in case of no pulse or only 1 pulse exists, 0.
    for j = 1:360/n
        tempopos = zeros( 1 , n); % temporary matrices to hold value of i so that the corresponding column in phi matrix can be easily found
        % by using phi(k,tempopos(1,x)).
        temponeg = zeros (1 , n);
        temposort = zeros (1, (3*n));% to sort the column values stored in tempopos and temponeg
        temposort2 = zeros (1, (3*n));
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
        pospulse = posi; % the total number of positive pulses in the window. this is for displaying.
        negpulse = nega;
        numpuls = posi + nega ;
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
                    deltat1 = ((-1)*deltat1)
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
                        deltat2 = ((-1)*deltat2)
                    end
                end
            else % when only 1 positive pulse is available
                disp (' There are no homogeneous pulses!! ' )
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
                    deltat1 = ((-1)*deltat1)
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
                        deltat2 = ((-1)*deltat2)
                    end
                end
            else
                disp (' There are no homogeneous pulses!! ' )
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
                    deltat1 = ((-1)*deltat1)
                end
                disp (' There are no homogeneous pulses!! ')
                remarks(k,j)=1;
                break
            else
                temposort = horzcat ( tempopos , temponeg)
                temposort2 = sort(temposort)
                temposort2(temposort2==0) = []; % to cut all the 0s which exist in the beg.
                %since the column values are arranged in ass. order thus leaving on,y the non-zero column values
                t1 = phi(k, temposort2(1,1));
                t2 = phi(k, temposort2(1,2));
                t3 = phi(k, temposort2(1,3));
                q1 = q(k, temposort2(1,1));
                q2 = q(k, temposort2(1,2));
                deltat1 = t2 - t1;
                q3 = q2 - q1;
                if ( q3 < 0)
                    deltat1 = ((-1)*deltat1)
                end
                deltat2 = t3 - t1;
                q3 = q(k, temposort2(1,3));
                q2 = q3 - q1;
                if ( q2 < 0)
                    deltat2 = ((-1)*deltat2)
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
            %disp (' There are no pulses!!! ')
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

