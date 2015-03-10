for i=1:ex
    for j=1:Nsamp
        eq_dist(j)=sqrt((w(i,:)-testdata(j,:))*((w(i,:)-testdata(j,:))'));
    end
    [temp,j1]=min(eq_dist);
    g(i,:)=testdata(j1,:);
    % store the sample which is at a minimum distance
end

%%%%%------calc observ-%%%%%%%

counter1=0;
%
% for j=1:ex
%     k=1;
%     for i=1:3:216
%         PHInew(j,k)=g(j,i);
%         Qnew(j,k)=g(j,i+1);
%         if Qnew(j,k)~=0
%             counter1=counter1+1;
%         end
%         Nnew(j,k)=g(j,i+2);
%         k=k+1;
%     end
% end

for j=1:ex
    k=1;
    for i=1:2:143
        PHInew(j,k)=g(j,i);
        Qnew(j,k)=g(j,i+1);
        if Qnew(j,k)~=0
            counter1=counter1+1;
        end
        %Nnew(j,k)=g(j,i+2);
        k=k+1;
    end
end


zz=k-1;
%clear Qpmod; % doubtful


for i=1:ex
    l=1;d=1;
    for j=1:144              % for j=1:216 if we use N also
        g1(i,l,d)=g(i,j);
        d=d+1;
        if d==3          % d==4 if we use N also
            d=1;
            l=l+1;
        end
    end
end

