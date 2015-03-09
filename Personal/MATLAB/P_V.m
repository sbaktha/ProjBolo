%A matlab function is written. The name of the function is given as 'P_V'
%because the file name and the function name must be the same. 'status'
%matrix is returned after the function body is executed.

function [ status ]= P_V()
global sum;global N;global i;global p_v;global j;global E_Pi;global E_A;global status;global tt;
sum = 0;
disp(' The probability of the node being visited during the training phase');
disp(N);
 for i=1:N 
     if(i==1)
        p_v(i)=E_Pi(i); 
     else 
        sum=0;
        for j=1:(i-1)
            sum= sum + p_v(j)*(E_A(j,i)/(1-E_A(j,j)) );      
        end
     end
                  p_v(i)= sum + (E_Pi(i));                      
 end
 
 tt=1;  
 for i=1:N
    if(p_v(i)*100 >= 40.0) 
       status(tt)=i;
       tt= tt +1;
    else
        status(i)=0;
    end
 end
end

          