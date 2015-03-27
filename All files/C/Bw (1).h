
void Bw_algo()
{      int kk=0;
        double sum2=0.0;

       // Calculation of ZI values.
       for(t=0;t<T-1;t++)
     {
          for(i=0;i<N;i++)
          {
               for(j=0;j<N;j++)
                {
                        nu=Alpha[t][i]*b[j][Ob[t+1]]*Beta[t+1][j]*a[i][j];
                        sum=0.0;
                        for( m=0;m<N;m++)
                        {
                           for( n=0;n<N;n++)
                           {
                              sum+=Alpha[t][m] *a[m][n]*b[n][Ob[t+1]] *Beta[t+1][n];
                           }  
                        }
                        ZI[t][i][j] =nu/sum;
                }      
            }
     }
     printf("\n ZI matrix is\n ");
     for(t=0;t<T-1;t++)
     {
          for(i=0;i<N;i++)
          {
               for(j=0;j<N;j++)
                {
                        printf("%lf \t",ZI[t][i][j]);
                } 
                printf("\n");                           
          }
          printf("\n");
       }
        // Gamma computation ..
        printf("GAMMA is");
        for(t=0;t<T;t++)
        {
           for(i=0;i<N;i++)
           {  sum=0.0;
              for(j=0;j<N;j++)
                 sum=sum+ZI[t][i][j];
              printf("%lf\t",Gamma[t][i]=sum);
           }     
           printf("\n");
        }   
        // Expected number of transistions from state i 
        for(i=0;i<N;i++)
        {   sum=0.0;
            for(t=0;t<T-1;t++)
               sum+=Gamma[t][i];
            E_T[i]=sum;
        }
        // Expected number of transitions from node i to node j 
      /*  for(i=0;i<N;i++)
        {
           for(j=0;j<N;j++)
           {   sum=0.0;
               for(t=0;t<T-1;t++)
                  sum=sum+ZI[t][i][j];
               E_I_J[i]=sum;                     // may be a mistake..............
           }       
        }   */
          // Computing estimated values for Pi ,A and B.
          printf("matrix pi is \n");
          for(i=0;i<N;i++)   //  for Pi
                printf("%lf \t ",E_Pi[i]=Gamma[0][i]);
           printf("\n");      
           for(i=0;i<N;i++)
           {
                 for(j=0;j<N;j++)
                 {
                         sum=0.0;nu=0.0;
                         for(t=0;t<T-1;t++)
                         {
                              sum=sum+ZI[t][i][j];
                              nu=nu+Gamma[t][i];
                         }   
                         E_A[i][j]= sum / nu  ;  
                 }        
           }      
          printf("The estimated transistion matrix is \n");
          for(i=0;i<N;i++)
           {
                 for(j=0;j<N;j++)
                 {
                     printf("\t %lf",E_A[i][j]);
                 }   
                 printf("\n");
           }      
          // Computing the matrix B.....
          printf("Matrix B is \n");
          for(j=0;j<N;j++)     // number of states 
          {   sum2=0.0; 
              for(kk=0;kk<K;kk++)
                  sum1[kk]=0.0;
              for(t=0;t<T;t++)   // to traverse the observation sequence...
              {   
                   for(kk=0;kk<K;kk++)
                   {
                      if(Ob[t] == kk)  // here one for loop will come
                      {    sum2+=Gamma[t][j]; // overall sum ..........
                           sum1[kk]+=Gamma[t][j];
                           break;
                      }
                   }
              }                                //here two more for loops will come
              
              for(kk=0;kk<K;kk++)
                  E_B[j][kk]= sum1[kk] /(sum2);
              for(kk=0;kk<K;kk++)    
                  printf("\t%lf \t",E_B[j][kk]);
               printf("\n");
          } 
     
}
