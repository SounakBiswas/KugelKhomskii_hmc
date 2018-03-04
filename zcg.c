#include "global.h"
//Solves Ax=b;
void zcg_init(double nd,double *x, double *b, int * rci, int *ipar,double *dpar,double *tmp){
  ipar[0]=nd;
  ipar[2]=1;
  ipar[3]=0;
  ipar[4]=nd;
  ipar[7]=0;
  ipar[8]=0;
  ipar[9]=1;
  ipar[10]=0;

  int i;
  for(i=0;i<4*nd;i++)
    tmp.real=tmp.imag=0;


}
void zcg(double nd,double *x, double *b, int * rci, int *ipar,double *dpar,double *tmp){
  niter=0;
  dcomplex cone,cmone,tvar;
  cone.real=1;cone.imag=-1;
  cmone.real=-1; cmone.imag=-1;
  tvar.imag=0;

  if(ipar[2]==1){
    zcopy(&nd,sol,&ione,tmp,&ione);
    zcopy(&nd, rhs, &ione, tmp+2*nd, &ione);
    *rci=1;
    ipar[2]=2;
    return;
    //tmp_nd=Ax0
  }
  else if(ipar[2]==2){
    zaxpy(&nd, &cmone, tmp+nd, &ione, tmp+2*nd, &ione);
    zcopy(&nd, tmp+2*nd, &ione, tmp, &ione);
    dpar[4]=zdotc(&nd, tmp+2*nd, &ione, tmp+2*nd, &ione );
    ipar[2]=3;
    *rci=1;
    return;

  }

  else if(ipar[2]==3){
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
    ipar[3]++;
    dpar[4]=dpar[5];
    *rci=2;
    return;
  }



}

