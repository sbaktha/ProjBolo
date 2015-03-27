%A matlab function is written. The name of the function is given as 'forw'
%because the file name and the function name must be the same. Alpha is
%returned after the function body is executed.
%The function call will be of the form : Alpha = forw;

function Alpha = forw()
global N;global Pi;global Ob;global b;global T;global a;
for i=1:N
        Alpha(1,i)=Pi(i) * b(i,Ob(1)); %Alpha(0,i),  Ob(0)
end
   for t=1:T-1
     for j=1:N
         sum=0;
         for i=1:N
           sum= sum+ Alpha(t,i)*a(i,j);
         end
          Alpha(t+1,j)=sum * b(j,Ob(t+1));
      end
   end
end