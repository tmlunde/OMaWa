#ifndef max
#define max(a,b) (((a) > (b)) ? (a) : (b))
#define min(a,b) (((a) < (b)) ? (a) : (b))
#endif

#include <stdlib.h>
#include <math.h>
#include <R.h>


void arabMort(double *tMax, int *nR, double *cCof, 
              double *d, int *mR, double *m1, double *m2, double *m3, double *m4, double *m5, 
	      double *m6, double *m7, double *m8, double *m9, double *mrt, double *bovine, double *human, double *HBI)
 	   
{
  int n;
  int m;
  n = *nR ;
  m = *mR ;





  /*local vars*/
  int i;
  int j;
  static double alpha;
  static double hbi;
  hbi = *HBI;
  
  for (i = 0; i <= (*nR - 1); ++i) {
   alpha = 0.1675256 + -0.0121402 * ((12 + pow((1+(tMax[i]+18)/11.10), (0.1686)) *
           (pow((1+(tMax[i]+18)/11.10), 1.991)  - (1+(tMax[i]+18)/11.10)*1.881 - 23))  - 
           250*(1*pow(100, (tMax[i]/3.))/4.641589e+26)) + exp(-(tMax[i] / 5.)) * 3. ;
   
   for (j = 0; j <= (*mR - 1); ++j) {
     mrt[j] = (1 + alpha*d[j])*exp(-alpha*d[j]) ;
   }
   
   m1[i] = (mrt[0] - mrt[1])/cCof[0] ;
   /*If arab does not find bloodmeal during day 2-4 mortality increases based on density of humans and bovine */
   m2[i] = min(hbi*((mrt[2] - mrt[4])/cCof[1] + (1-human[i])) + (1-hbi)*((mrt[2] - mrt[4])/cCof[1] + (1-bovine[i])), 1);
   m3[i] = (mrt[5] - mrt[8])/cCof[2] ;
   m4[i] = (mrt[9] - mrt[13])/cCof[3];
   m5[i] = (mrt[14] - mrt[19])/cCof[4];
   m6[i] = (mrt[20] - mrt[26])/cCof[5];
   m7[i] = (mrt[27] - mrt[34])/cCof[6];
   m8[i] = (mrt[35] - mrt[43])/cCof[7];
   m9[i] = (mrt[44] - mrt[199])/cCof[8];

  }
   
           
  
  
}






