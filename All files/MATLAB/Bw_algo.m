%A matlab function is written. The name of the function is given as 'Bw_algo'
%because the file name and the function name must be the same. 'ZI , Gamma, E_Pi , E_A' and
%'E_B' are returned after the function body is executed.

function [ZI , Gamma, E_Pi , E_A , E_B ]= Bw_algo()
global kk; global sum2; global t;global T;global ZI; global Gamma;
global E_Pi; global E_A; global E_B; global N; global m;global n;
global nu;global sum;global a;global b;global Ob;global Beta;
global E_T;global K;global sum1;global Alpha;global i;global j;
kk=0;
sum2=0;
%Calculation of ZI values
for t=1:T-1
    for i=1:N
        for j=1:N
            nu=Alpha(t,i)*b(j,Ob(t+1))*Beta(t+1,j)*a(i,j);
            sum=0;
            for m=1:N
                for n=1:N
                    sum = sum + (Alpha(t,m) *a(m,n) *b(n,Ob(t+1)) *Beta(t+1,n));
                end  
            end
                        ZI(t,i,j) = nu/sum;
        end      
    end
end
    
 %Gamma computation
    for t=1:T    
        for i=1:N
            sum=0;
            for j=1:N
                sum = sum+ZI(t,i,j);
            end
            Gamma(t,i)=sum;
        end     
    end
    
  %Expected number of transistions from state i
        for i=1:N
            sum=0;
            for t=1:T-1
               sum= sum + Gamma(t,i);
            end
            E_T(i)=sum;
        end
        
   %Expected number of transitions from node i to node j 
      for i=1:N
           for j=1:N
               sum=0;
               for t=1:T-1
                   sum= sum+ZI(t,i,j);
               end
               E_I_J(i)=sum;                     % may be a mistake (already mentioned in the C version)..............
           end       
      end
      
     %Computing estimated values for Pi ,A and B.
          for i=1:N   
              E_Pi(i)= Gamma(1,i);    % E_Pi(i)= Gamma(0,(i));
          end
          
          for i=1:N
              for j=1:N
                  sum=0;
                  nu=0;
                  for t=1:T-1
                      sum=sum+ZI(t,i,j);
                      nu=nu+Gamma(t,i);
                  end   
                         E_A(i,j) = (sum / nu)  ;  
              end        
          end
          
       %Computing the matrix B
            for j=1:N     % number of states 
                sum2=0; 
                for kk=1:K
                    sum1(kk)=0;
                end
                for t=1:T  %to traverse the observation sequence...
                    for kk=1:K
                        if(Ob(t) == kk)  % here one for loop will come
                           sum2 = sum2+ Gamma(t,j); % overall sum ..........
                           sum1(kk)= sum1(kk) + Gamma(t,j);
                           break;
                        end
                    end
                end
                for kk=1:K
                  E_B(j,kk) = (sum1(kk))/sum2;
                end
            end
end