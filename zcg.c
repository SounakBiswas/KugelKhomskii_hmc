#include "global.h"
//Solves Ax=b;
void zcg(double nd,double *x, double *b, int * rci, double *dpar,double *tmp){
  niter=0;
  dcomplex cone,cmone,tvar;
  cone.real=1;cone.imag=-1;
  cmone.real=-1; cmone.imag=-1;
  tvar.imag=0;


  zcopy(&nd,sol,&ione,tmp,&ione);
  //tmp_nd=Ax0
  zcopy(&nd, rhs, &ione, tmp+2*nd, &ione);
  zaxpy(&nd, &cmone, tmp+nd, &ione, tmp+2*nd, &ione);
  zcopy(&nd, tmp+2*nd, &ione, tmp, &ione);
  dpar[5]=zdotc(&nd, tmp+2*nd, &ione, tmp+2*nd, &ione );
  dpar[4]=dpar[5];
  *rci=1;
  return;


  //tmp_nd=A tmp
  dpar[6]=dpar[5]/zdotc(&nd, tmp+nd, &ione, tmp+2*nd, &ione );
  tvar.real=dpar[6];
  zaxpy(&nd, &tvar, tmp, &ione, sol, &ione );
  tvar.real=-dpar[6];
  zaxpy(&nd, &tvar, tmp+nd, &ione, tmp+2*nd, &ione );
  dpar[5]=zdotc(&nd, tmp+2*nd, &ione, tmp+2*nd, &ione );
  dpar[7]=dpar[5]/dpar[4];
  tvar.real=dpar[7];
  zdscal(&nd, &dpar[7], tmp, &ione);
  zaxpy(&nd, &cone, tmp+2*nd, &ione, tmp+nd, &ione );
  niter++;
  dpar[4]=dpar[5];
  *rci=2;
  return;



}

