void  forward()
{
      //Forward Calculation 
     for(i=0;i<N;i++)
        Alpha[0][i]=Pi[i] * b[i][Ob[0]];
     for(t=0;t<T-1;t++)
     {
        for(j=0;j<N;j++)
        {
           sum=0.0; 
           for(i=0;i<N;i++)
               sum=sum+ Alpha[t][i]*a[i][j];
           Alpha[t+1][j]=sum * b[j][Ob[t+1]];
        }   
     }   
     printf("Forward matrix is \n");// printing Forward matrix ....
     for(t=0;t<T;t++)
     {
         for(i=0;i<N;i++)
           printf("%lf\t",Alpha[t][i]);
         printf("\n");
     }
}
