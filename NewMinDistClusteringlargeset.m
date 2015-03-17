%----------------------------%
z=0;
ex = 18;
%testdata(1:Nsamp,:) = phiqdata(1:Nsamp,:);

PHIpmod=PHIp;

for i=1:Nsamp
    kk=1;
    for j=1:72
        testdata(i,kk)=PHIpmod(i,j);
        testdata(i,kk+1)=Qpmod(i,j);
        kk=kk+2;
    end
end

% testdata(1:Nsamp,:)=Qpmod(1:Nsamp,:);

%%%%modified here.. removed data, tot,testdata.
k=1;
for i=1:18
    w(i,:)=testdata(k,:);
    k=k+floor(Nsamp/18);
end

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

