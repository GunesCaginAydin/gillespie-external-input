#include <stdio.h>
#include <math.h>
#define __HEADERS__
#include "cowan.c"

int main(int argc,char **argv)
  {
   int N=100000;
   double ws=13.8;
   double w0=0.1;
   double h=1e-6;

   double alpha=0.1;
   double beta=1;
   double delta=0.1;		// step temporale in millisecondi

   int cmin=1000;
   int cmax=10000;

   cowan *X=new cowan(N,h,ws,w0,alpha,beta);

   int imax=40;
   int *del=(int*)calloc(imax,sizeof(int))-1;
   del[1]=1;
   for (int i=2;i<=imax;i++) del[i]=(int)ceil(del[i-1]*1.122);
   int tmax=del[imax];

   double B=0;
   double B2=0;
   double *C=(double*)calloc(imax,sizeof(double))-1;

   for (int c=1;c<=cmax;c++)
     {
      double f0=X->ftmp;
      int i=1;
      for (int t=1;t<=tmax;t++)
        {
         X->dinamica(delta);
         if (c>cmin)
           {
            double ftmp=X->ftmp;
            B+=ftmp;
            B2+=ftmp*ftmp;
            if (t==del[i])
              {
               C[i]+=f0*ftmp;
               i++;
              }
           }
        }
      if (c%1000==0) printf("c=%i/%i\n",c,cmax);
     }

   double bm=B/(tmax*(cmax-cmin));
   double b2m=B2/(tmax*(cmax-cmin));
   for (int i=1;i<=imax;i++) printf("%12i %12g\n",del[i],(C[i]/(cmax-cmin)-bm*bm)/(b2m-bm*bm));

   return 0;
  }
