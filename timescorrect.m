clear all;
for aaaa=1:9
    clearvars -except aaaa
    fileno=strcat('00',num2str(aaaa));
    type='Surface';
    filename1=strcat('J:\Datasets\',type,'\',type,'File',fileno);
    phase=textread(strcat(filename1,'.pdb.P.txt'));
    times=textread(strcat(filename1,'.pdb.Ti.txt'));
    offset=phase(1,1)*0.02/360;
    dlmwrite(strcat(filename1,'.pdb.Tioffset.txt'),(times+offset),'newline','pc')
end
