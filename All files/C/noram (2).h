void noramal()
{    double sum2=0.0,sum3=0.0;
     int pp=0,pp1=0;
     for(i=0;i<N;i++)   // for each row ....
     {
         if(i==status[pp])  // row is not normalised  
         {     pp=pp+1;
               for(j=0;j<N;j++)
                 N_E_A[i][j]=E_A[i][j];     // row copied suc
               continue ;
         }
         else              // row is  normalised 
         {
             sum3=0.0;sum2=0.0;pp1=0;
             for(j=0;j<N;j++)      // for each column 
             {
                 if(j==status[pp1]) // get the sum of edge weights keeping as it is ..
                  { pp1=pp1+1;
                    sum3+=a[i][j];
                  }
                 else 
                    sum2+=E_A[i][j];  //get the sum of edge weights to be normalised w.r.t sum3 ...
             }
             pp1=0;
             for(j=0;j<N;j++)
             {
                 if(j!=status[pp1])
                 {
                       N_E_A[i][j]=(1-sum3)*(E_A[i][j]/sum2);
                 }
                 else 
                 {
                         pp1=pp1+1;
                         N_E_A[i][j]=a[i][j];
                 }  
             }
         }
     }
     printf("\nBEFORE NORMALISATION\n");
     for(i=0;i<N;i++)
     {
        for(j=0;j<N;j++)
            printf("%lf\t",E_A[i][j]);
        printf("\n");
     }
     printf("\nAFTER NORMALISATION\n");
     for(i=0;i<N;i++)
     {
        for(j=0;j<N;j++)
            printf("%lf\t",N_E_A[i][j]);
        printf("\n");
     }
     
                    



}
                             
                 
             
