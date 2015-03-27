void P_V()      
{
          int i,j;
          double sum=0;
          // skip train .....
          
          printf("The probability of the node being \n\tvisited during the training phase %d\n",N);
         
          for(i=0;i<N;i++)  // loop for each node 
          {   
                  if(i==0)
                  { p_v[i]=E_Pi[i]; 
                    continue;                   
                  }
                  else 
                  {    sum=0;
                       printf("\n");
                       for(j=0;j<=i-1;j++)
                          sum=sum+ p_v[j]*(E_A[j][i]/(1-E_A[j][j]) );      
                  }
                  p_v[i]= sum + (E_Pi[i]);                      
                  
           }
           printf("\n");
           for(i=0;i<N;i++)
             printf("%lf\t",p_v[i]*100);
           tt=0;
           for(i=0;i<N;i++)
           {
              if(p_v[i]*100 >=40.0) 
              {
                status[tt]=i;
                tt++;
              }  
           }    
              //else 
               // status[i]=0;
          printf("\nnstatus matrix is \n");
          for(i=0;i<tt;i++)
             printf("\t%d",status[i]);
          
             
}
