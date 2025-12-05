#ifndef __COWAN_H__
#define __COWAN_H__

class cowan
  {
   public:
   cowan(long Ne,long Ni,double ws,double w0,double h,double alpha,double beta);
   double fsigma(void);
   long dinamica(double step);

   long Ne;
   long Ni;
   double h;
   double w0;
   double ws;
   double wa;
   double wee;
   double wei;
   double wie;
   double wii;
   double alpha;
   double beta;
   long k;		// neuroni eccitatori attivi
   long l;		// neuroni inibitori attivi
  };

#ifndef __HEADERS__
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// find fixed point of the dynamics
double cowan::fsigma(void)
  {
   double s1=0;
   double s2=1;
   while (fabs(s2-s1)>1e-12)
     {
      double sm=(s1+s2)/2;
      double s=w0*sm+h;
      double rate=(s<=0?0:beta*tanh(s));
      double ftmp=(1-sm)*rate-alpha*sm;
      if (ftmp<0) s2=sm; else s1=sm;
     }
   return (s1+s2)/2;
  }

cowan::cowan(long _Ne,long _Ni,double _ws,double _w0,double _h,double _alpha,double _beta)
  {
   Ne=_Ne;	// copy parameters
   Ni=_Ni;
   h=_h;
   ws=_ws;
   w0=_w0;
   alpha=_alpha;
   beta=_beta;
   double S0=fsigma();
   k=(long)rint(Ne*S0);		// start the dynamics at the fixed point
   l=(long)rint(Ni*S0);
   wee=(ws+w0)/2;
   wie=(ws+w0)/2;
   wei=(ws-w0)/2;
   wii=(ws-w0)/2;
  }

// Gillespie dynamics
long cowan::dinamica(double step)
  {
   long nhit=0;
   double t=0;
   while (t<step)
     {
      double x=k/(double)Ne;
      double y=l/(double)Ni;
      double se=wee*x-wei*y+h;
      double si=wie*x-wii*y+h;
      double fe=(se<=0?0:beta*tanh(se));		// activation rate excitatory neurons
      double fi=(si<=0?0:beta*tanh(si));		// 	"	"  inhibitory	"
      double rae=(Ne-k)*fe;				// total activ. rate excitatory
      double rai=(Ni-l)*fi;				// total activ. rate inhibitory
      double rse=k*alpha;				// total deactiv. rate excitatory
      double rsi=l*alpha;				// total deactiv. rate inhibitory
      double rtot=rae+rai+rse+rsi;
      double dt=-log((1.+random())/(2.+RAND_MAX))/rtot;		// temporal step before next event
      if (t+dt>step) break;
      t+=dt;
      double xtmp=rtot*(1.+random())/(2.+RAND_MAX);
      if (xtmp<rae) {k++; nhit++;}
      else if (xtmp<rae+rai) {l++; nhit++;}
      else if (xtmp<rae+rai+rse) k--;
      else l--;
     }
   return nhit;						// number of spikes in the interval
  }
#endif
#endif
