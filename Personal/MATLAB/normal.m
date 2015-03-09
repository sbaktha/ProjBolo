%A matlab function is written. The name of the function is given as 'noramal'
%because the file name and the function name must be the same. 'N_E_A' and
%'E_A' are returned after the function body is executed.

function [ N_E_A , E_A ]= normal()
global sum2;global sum3;global pp;global pp1;global i;global N;global N_E_A;global E_A;global j;global status;global a;
sum2 = 0;
sum3 = 0;
pp = 0;
pp1 = 0;
for i=1:N
    if(i==status(pp+1)) % status(pp)
        pp=pp+1;
         for j=1:N
         N_E_A(i,j)=E_A(i,j);
         end
    else
        sum3=0;
        sum2=0;
        pp1=0;
             for j=1:N       
                 if(j==status(pp1+1)) %status(pp1)
                    pp1=pp1+1;
                    sum3= sum3 + a(i,j);
                 else 
                    sum2= sum2 + E_A(i,j);
                 end
             end
             
             pp1=0;
             for j=1:N
                 if(j~=status(pp1+1))  %status(pp1)
                   N_E_A(i,j)=(1-sum3)*(E_A(i,j)/sum2);
                 else
                   pp1=pp1+1;
                   N_E_A(i,j)=a(i,j);
                 end
             end
    end
end
end
