b1=zeros(5,72);
for l=1:o
b1(:,:)=B(l,:,:);
end



o=18;
for i=1:o
l=1;
for j=1:72
for k=1:st
b2(i,l)=B(i,k,j);
l=l+1;
end
end
end


for i=1:o
for j=1:72
b3(i,j)=B(i,statexy(i,j),j);
end
end