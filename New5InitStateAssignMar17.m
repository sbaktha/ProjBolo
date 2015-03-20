
%%check if i should use chi-means distance and not normal euclidean as the
%%effect of PHI is more in the state assignment.
st=4;

for k=1:ex
    clear g;
    g(:,:)=g1(k,:,:);
    for i=1:st
        c(i,:)=g(i,:);
    end
    cl=1;
    coun=1;
    while cl==1
        % cl=1;
        %coun
        for i=1:st
            for j=1:wind
                dm(i,j)=sqrt((c(i,:)-g(j,:))*((c(i,:)-g(j,:))'));
            end
        end
        sm=zeros(st,wind);
        [temp,temp1]=min(dm);
        for i=1:wind
            sm(temp1(i),i)=1;
        end
        tsm=sum(sm,2);
        for i=1:st
            if tsm(i,1)>1
                temp2=0;
                for j=1:wind
                    if sm(i,j)==1
                        temp2=temp2+g(j,:);
                    end
                end
                c(i,:)=temp2/tsm(i,1);
            end
        end
        
        for i=1:st
            for j=1:wind
                dm(i,j)=sqrt((c(i,:)-g(j,:))*((c(i,:)-g(j,:))'));
            end
        end
        
        sm1=zeros(st,wind);
        [temp,temp1]=min(dm);
        for i=1:wind
            sm1(temp1(i),i)=1;
        end
        % coun=coun+1;
        if sm==sm1
            cl=0;
        end
        
    end
    sumstate(:,k)=sum(sm1,2);
    ss(k,:)=temp1(1,:);
end

S=ss;


%%%%%Old
% wind=72;
% ex=18;
% for k=1:ex
%     clear g;
%     g(:,:)=g1(k,:,:);
%     st=4;
%     for i=1:st
%         c(i,:)=g(i,:);
%     end
%     cl=1;
%     while cl==1
%         cl=1;
%         for i=1:st
%             for j=1:wind
%                 dm(i,j)=sqrt((c(i,:)-g(j,:))*((c(i,:)-g(j,:))'));
%             end
%         end
%         sm=zeros(st,wind);
%         [temp,temp1]=min(dm);
%         for i=1:wind
%             sm(temp1(i),i)=1;
%         end
%         tsm=sum(sm,2);
%         for i=1:st
%             if tsm(i,1)>1
%                 temp2=0;
%                 for j=1:wind
%                     if sm(i,j)==1
%                         temp2=temp2+g(j,:);
%                     end
%                 end
%                 c(i,:)=temp2/tsm(i,1);
%             end
%         end
%         for i=1:st
%             for j=1:wind
%                 dm(i,j)=sqrt((c(i,:)-g(j,:))*((c(i,:)-g(j,:))'));
%             end
%         end
%         sm1=zeros(st,wind);
%         [temp,temp1]=min(dm);
%         for i=1:wind
%             sm1(temp1(i),i)=1;
%         end
%         if sm==sm1
%             cl=0;
%         end
%     end
%     sumstate(:,k)=sum(sm1,2);
%     ss(k,:)=temp1(1,:);
% end
%
% S=ss;

