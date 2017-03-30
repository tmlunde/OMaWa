#ifndef max
#define max(a,b) (((a) > (b)) ? (a) : (b))
#define min(a,b) (((a) < (b)) ? (a) : (b))
#endif

#include <stdlib.h>
#include <math.h>
#include <R.h>


void gambMort(double *temp, double *RH, int *nR, double *cCof, 
              double *d, int *mR, double *m1, double *m2, double *m3, double *m4, double *m5, 
	      double *m6, double *m7, double *m8, double *m9, double *mrt, 
	      double *bovine, double *human, double *HBI, double *f, int *kR, double *myK, double *mSiz)
 	   
{
  int n;
  int m;
  n = *nR ;
  m = *mR ;
  /*local vars*/
  int i;
  int j;
  int k;
  static double alpha;
  static double nlsEstCoef;
  static double cVal;
  static double cbVal;
  static double hbi;
  static double hM;
  hbi = *HBI;
  cVal = 21;
  cbVal = -1;
  
  for (i = 0; i <= (*nR - 1); ++i) {
   hM = 2.1731 - 0.3846 * mSiz[i]  ;
   nlsEstCoef = 6.48007 + 0.69570*(1-exp(-0.06*RH[i])) ;
   alpha = (exp(10 +  pow((1+(temp[i]-cbVal)/cVal), 0.6666667) *
        ((1+(temp[i]-cbVal)/cVal)*(1+(temp[i]-cbVal)/cVal) - (1+(temp[i]-cbVal)/cVal)*2. - nlsEstCoef)))*hM;
   
   /*Rprintf("Error: %f \n", alpha) ;*/
   
   for (j = 0; j <= (*mR - 1); ++j) {
     mrt[j] = 0 ;
     for (k = 0; k <= (*kR - 1); ++k) {
       mrt[j] += (pow(alpha*d[j], myK[k])/f[k])*exp(-alpha*d[j]) ;
     }
   }
   
   
   m1[i] = -log(mrt[1]/mrt[0])/cCof[0] ;
   /*If arab does not find bloodmeal during day 2-4 mortality increases based on density of humans and bovine */
   m2[i] = min(hbi*(-log(mrt[4]/mrt[2])/cCof[1] + (1-human[i])) + (1-hbi)*(-log(mrt[4]/mrt[2])/cCof[1] + (1-bovine[i])), 1);
   m3[i] = -log(mrt[8]/mrt[5])/cCof[2] ;
   m4[i] = -log(mrt[13]/mrt[9])/cCof[3];
   m5[i] = -log(mrt[19]/mrt[14])/cCof[4];
   m6[i] = -log(mrt[26]/mrt[20])/cCof[5];
   m7[i] = -log(mrt[34]/mrt[27])/cCof[6];
   m8[i] = -log(mrt[43]/mrt[35])/cCof[7];
   m9[i] = -log(mrt[199]/mrt[44])/cCof[8];

  }
  
  
  
   
           
  
  
}
