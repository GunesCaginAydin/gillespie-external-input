#ifndef __COWAN_H__
#define __COWAN_H__

class cowan
  {
   public:
   cowan(long N,double h,double ws,double w0,double alpha,double beta);
   double fsigma(void);
   void start(long *k0,long *l0);
   double dinamica(double step);

   long N;
   double h;
   double w0;
   double ws;
   double alpha;
   double beta;
   long k;		// neuroni eccitatori attivi
   long l;		// neuroni inibitori attivi
   double ftmp;
  };

#ifndef __HEADERS__
#include <stdio.h>
#include <math.h>

cowan::cowan(long _N,double _h,double _ws,double _w0,double _alpha,double _beta)
  {
   N=_N;
   h=_h;
   ws=_ws;		// w_E+w_I
   w0=_w0;		// W_E-w_I
   alpha=_alpha;
   beta=_beta;
   double S0=fsigma();
   k=(long)rint(S0*N);
   l=(long)rint(S0*N);
   double sigma=(k+l)/(2.*N);                           // frazione di neuroni attivi
   double delta=(k-l)/(2.*N);                           // frazioni di neuroni eccitatori - neuroni inibitori
   double s=w0*sigma+ws*delta+h;                                // campo locale sul neurone
   ftmp=2*(1-sigma)*(s<=0?0:beta*tanh(s));              // firing rate istantaneo
  }

double cowan::fsigma(void)
  {
   double s1=0;
   double s2=1;
   while (fabs(s2-s1)>1e-12)
     {
      double sm=(s1+s2)/2;
      double s=w0*sm+h;
      double rate=(s<=0?0:beta*tanh(s));
      double fm=(1-sm)*rate-alpha*sm;
      if (fm<0) s2=sm; else s1=sm;
     }
   return (s1+s2)/2;
  }

double cowan::dinamica(double step)
  {
   int nhit=0;
   double t=0;
   while (t<step)
     {
      double sigma=(k+l)/(2.*N);
      double delta=(k-l)/(2.*N);
      double s=w0*sigma+ws*delta+h;
      double rate=(s<=0?0:beta*tanh(s));		// rate di accensione dei singoli neuroni
      double rae=(N-k)*rate;				// rate totale accensione eccitatori
      double rai=(N-l)*rate;				// rate totale accensione inibitori
      double rse=k*alpha;				// rate totale spegnimento eccitatori
      double rsi=l*alpha;				// rate totale spegnimento inibitori
      double rtot=rae+rai+rse+rsi;
      double dt=-log((1.+random())/(2.+RAND_MAX))/rtot;		// step temporale prima del prossimo evento
      if (t+dt>step) break;
      t+=dt;
      double xtmp=rtot*(1.+random())/(2.+RAND_MAX);
      if (xtmp<rae) {k++; nhit++;}
      else if (xtmp<rae+rai) {l++; nhit++;}
      else if (xtmp<rae+rai+rse) k--;
      else l--;
     }
   //return nhit;						// numero di spike nell'intervallo
   double sigma=(k+l)/(2.*N);				// frazione di neuroni attivi
   double delta=(k-l)/(2.*N);				// frazioni di neuroni eccitatori - neuroni inibitori
   double s=w0*sigma+ws*delta+h;				// campo locale sul neurone
   ftmp=2*(1-sigma)*(s<=0?0:beta*tanh(s));		// firing rate istantaneo
   return ftmp;
  }
#endif
#endif
