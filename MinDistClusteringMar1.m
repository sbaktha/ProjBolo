%----------------------------%
z=0;
ex = 18;
%testdata(1:Nsamp,:) = phiqdata(1:Nsamp,:);

PHIpmod=PHIp;

for i=1:30
    kk=1;
    for j=1:72
        testdata(i,kk)=PHIpmod(i,j);
        testdata(i,kk+1)=Qpmod(i,j);
        kk=kk+2;
    end
end

% testdata(1:Nsamp,:)=Qpmod(1:Nsamp,:);

%%%%modified here.. removed data, tot,testdata.
w(1,:) = testdata(1,:);
w(2,:)= testdata(3,:);    % i've used only Q data as
w(3,:) = testdata(5,:);   % instructed by prof.
w(4,:)= testdata(6,:);
w(5,:) = testdata(8,:);
w(6,:) = testdata(10,:);
w(7,:)= testdata(11,:);
w(8,:)= testdata(13,:);
w(9,:) = testdata(15,:);
w(10,:) = testdata(16,:);
w(11,:) = testdata(18,:);
w(12,:) = testdata(20,:);
w(13,:) = testdata(21,:);
w(14,:) = testdata(23,:);
w(15,:) = testdata(25,:);
w(16,:) = testdata(26,:);
w(17,:) = testdata(28,:);
w(18,:) = testdata(30,:);

% Columns 1 through 15
%1     4     2     4     3     4     7     5     5     6     7     1     8     7     9
%Columns 16 through 30
% 10     2    11    12    12    13     7    14     7    15    16    14    17     6    18

%
% w(1,:) = Qpmod(1,:);
% w(2,:)= Qpmod(3,:);
% w(3,:) = Qpmod(5,:);
% w(4,:)= Qpmod(6,:);
% w(5,:) = Qpmod(8,:);
% w(6,:) = Qpmod(10,:);
% w(7,:)= Qpmod(11,:);
% w(8,:)= Qpmod(13,:);
% w(9,:) = Qpmod(15,:);
% w(10,:) = Qpmod(16,:);
% w(11,:) = Qpmod(18,:);
% w(12,:) = Qpmod(20,:);
% w(13,:) = Qpmod(21,:);
% w(14,:) = Qpmod(23,:);
% w(15,:) = Qpmod(25,:);
% w(16,:) = Qpmod(26,:);
% w(17,:) = Qpmod(28,:);
% w(18,:) = Qpmod(30,:);


w_o = zeros(size(w)); % create a zero matrix of same dimentions of weight
dif = w - w_o; % initialise difference.
count = 1;                   % initialise iteration count.
while sum(dif(:)) ~= 0 && count ~= 500
    sum(dif(:));
    w_o = w; % remember the weights of previous iterations.
    for ii = 1:Nsamp
        for jj = 1:ex
            eq_dist(jj) = ((testdata(ii,:)-w(jj,:)) * ((testdata(ii,:)-w(jj,:))')) ; % equiledian distance
        end
        [temp,near_class(ii)] = min(eq_dist);  % find the cluster which is in minimum distance from the training exempler.
    end
    for ii = 1:ex
        [a,b] = find(near_class == ii);
        temp_sum = 0;
        for jj = b
            temp_sum = temp_sum + testdata(jj,:);
        end
        if sum(a) == 0
            count;
            ii;
        end
        w(ii,:) = temp_sum / sum(a);
    end
    dif = abs(w - w_o);
    count = count+1;
    %------min dist----------%
end

