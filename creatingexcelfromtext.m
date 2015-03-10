clear all;
for i=1:9
    clearvars -except i;
    fileno=strcat('00',num2str(i));
    type='Corona';
    filename1=strcat('J:\Datasets\',type,'\',type,'File',fileno);
    ampli=textread(strcat(filename1,'.pdb.A.txt'));
    phase=textread(strcat(filename1,'.pdb.P.txt'));
    times=textread(strcat(filename1,'.pdb.Ti.txt'));
    %times=textread('J:\Universita di Bologna\Papers for modification in code\venkatesh sir 4 feb files\data of epoxy double v\epoxy void_ac test.024.pd2.Ti.txt');
    var(:,1)=phase;
    var(:,2)=ampli;
    var(:,3)=times;
    xlswrite(strcat(filename1,'-rawPQTIdata.xlsx'),var);
end

