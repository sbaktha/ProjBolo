void backward()
{
     // Back ward calculation .....
     for(i=0;i<N;i++)
        Beta[T-1][i]=1;
     for(t=T-2;t>=0;t--)
     {
          for(i=0;i<N;i++)
          {
              sum=0.0;
              for(j=0;j<N;j++)
                 sum=sum+a[i][j]*Beta[t+1][j]*b[j][Ob[t+1]];
              Beta[t][i]=sum;
          }   
     }    
     printf("Backward matrix is \n");                //Printing Backward matrix ...  
     for(t=0;t<T;t++)
     {
         for(i=0;i<N;i++)
           printf("%lf\t",Beta[t][i]);
         printf("\n");
     }    
}
