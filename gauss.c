#include "global.h"
#include <stdlib.h> 
#include <stdio.h> 
#include <math.h>
#include "mt19937ar.h"
// Generates Gaussian Random numbers
// Based on Marsaglia algorithm on wikipedia :https://en.wikipedia.org/wiki/Marsaglia_polar_method 

double genrand_gauss(double mean,double var){
  static double z1,z2;
  static int call=0;
  double u1,u2;
  double w;
  call= (call+1)%2;
 

  if(call) {
    while(1) {
      u1=-1+2*genrand_real2();
      u2=-1+2*genrand_real2();
      w=u1*u1+u2*u2;
      if(w<1)
        break;
    }
    z1=mean+var*u1*sqrt(-2.0*log(w)/w);
    z2=mean+var*u2*sqrt(-2.0*log(w)/w);
    return z1;
  }
  else
    return z2;

}
void rand_norm_vec(double *vec,int dim,double mean,double var){
  int i;
  for(i=0;i<dim;i++){
    vec[i]=genrand_gauss(mean,var);
  }
}
void rand_norm_cvec(dcomplex *vec,int dim,double mean,double var){
  int i;
  for(i=0;i<dim;i++){
    vec[i].real=genrand_gauss(mean,var);
    vec[i].imag=genrand_gauss(mean,var);
  }
}
