/*

The code writes on stdout the duration and size of the avalanches.
Avalanches are defined as a consecutive sequence of time bins where at least one spike is observed in the system.
The first cmin avalanches are skipped for thermalization, and the program stops when cmax avalanches have been counted.

*/
#include <stdio.h>
#include <math.h>
#define __HEADERS__
#include "cowan.c"

int main(int argc,char **argv)
  {
   double w0=0.1;			// difference of excitatory and inhibitory connections
   double h=1e-6;			// external input
   long N=1000000;			// number of excitatory neurons
   double delta=0.64;                   // time bin

   // parameters to be held constant
   double ws=13.8;			// sum of excitatory and inhibitory connections
   double alpha=0.1;			// deactovation rate
   double beta=1;			// activation rate coefficient of hyerbolic tangent
   double xin=0.5;			// fraction of inhibitory neurons
   long Ne=N;					// number of excitatory neurons
   long Ni=(int)rint(xin*N/(1-xin));		// number of inhibitory neurons
   cowan *X=new cowan(Ne,Ni,ws,w0,h,alpha,beta);

   long stmp=0;
   int n=0;
   long count=0;			// counts avalanches
   long cmin=10000;			// min number of avalanches
   long cmax=1000000;			// max number of avalanches
   while (1)
     {
      long nhit=X->dinamica(delta);
      if (nhit==0)
        {
         if (n>0)
           {
            // end of an avalanche
            count++;
            if (count>cmin) printf("%12i %12li\n",n,stmp);	// print size and duration of the avalanche
            if (count>=cmax) break;
            stmp=0;
            n=0;
           }
        }
      else 
        {
         stmp+=nhit;
         n++;
        }
     }

   return 0;
  }
