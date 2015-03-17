%plotting the Clusters
typepl=['m*';'c*';'r*';'g*';'b*';'k*';'mo';'co';'ro';'go';'bo';'ko';'m+';'c+';'r+';'g+';'b+';'k+'];
for k=7:18
    for i=1:2998
        if near_class(i)==k
            scatter(PHIp(i,:),Qpmod(i,:),10,typepl(near_class(i),:));
            hold on
        end
    end
    hold off;
    saveas(gcf,strcat('cluster',num2str(k),'.jpeg'));
    close Figure 1;
end
