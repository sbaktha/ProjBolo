%A matlab function is written. The name of the function is given as 'backw'
%because the file name and the function name must be the same. Beta is
%returned after the function body is executed.
%The function call will be of the form : Beta = backw;
function Beta = backw()
global i;global N;global T;global Beta;global t;global j;global a;global b;global Ob;global sum;
   for i=1:N
       Beta(T,i)=1; 
   end                   %Beta=ones(T,N);
   for t=T-1:-1:1        %t-T-2:1
     for i=1:N
         sum=0;
         for j=1:N
           sum=sum+(a(i,j)*Beta(t+1,j)*b(j,Ob(t+1)));
         end
          Beta(t,i)=sum;
      end
   end
end