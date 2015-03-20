clear all
clc
for iiii=1:1
    clc;
    clearvars -except iiii ;
    fileno=strcat('00',num2str(iiii));
    type='Internal';
    New1StartingCode;
    UniformInitforLRmodel;
    HMMCalcMar15;
    
    filenamemodel=strcat('Lambda',type,'file',num2str(fileno),'.mat');
    save(filenamemodel,'Pi','N_E_A','N_E_B','New_A','New_B','mlop','mlogprob');
    %delay(5);
end


