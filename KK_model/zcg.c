#include "global.h"
#include "stdio.h"
//Solves Ax=b;
void zcg_init(int nd,dcomplex *x, dcomplex *b, int * rci, int *ipar,double *dpar,dcomplex *tmp){
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
    tmp[i].real=tmp[i].imag=0;


}
void zcg(int nd,dcomplex *sol, dcomplex *rhs, int * rci, int *ipar,double *dpar,dcomplex *tmp){
  dcomplex cone,cmone,tvar;
  cone.real=1;cone.imag=0;
  cmone.real=-1; cmone.imag=0;
  tvar.imag=0;
  //printf("ipar[2]=%d\n",ipar[2]);

  //Initialization
  if(ipar[2]==1){
    ipar[3]=0;
    //p_0=x_0
    zcopy(&nd,sol,&ione,tmp,&ione);

    //r_0=b
    zcopy(&nd, rhs, &ione, tmp+2*nd, &ione);
    //calcluate &tmp[nd] <----  A.p
    *rci=1;
    ipar[2]=2;
    return;
    //tmp_nd=Ax0
  }
  else if(ipar[2]==2){
    //r_0=b-Ax
    zaxpy(&nd, &cmone, tmp+nd, &ione, tmp+2*nd, &ione);

    //p_0=r_0
    zcopy(&nd, tmp+2*nd, &ione, tmp, &ione);

    //Calc <r0,r0>
    zdotc(&tvar, &nd, tmp+2*nd, &ione, tmp+2*nd, &ione );
    dpar[4]=tvar.real;
    //printf("init norm :%f\n",dpar[4]);
    ipar[2]=3;
    //put Ap into temp+nd
    *rci=1;
    return;

  }

  else if(ipar[2]==3){
    // \alpha = <r,r>/<p,Ap>
    zdotc(&tvar, &nd, tmp, &ione, tmp+nd, &ione );
    dpar[6]=dpar[4]/tvar.real;
  //  printf("rold=%f\talpha=%f\n",dpar[4],dpar[6]);
    tvar.real=dpar[6];
    // x= x + \alpha p
    zaxpy(&nd, &tvar, tmp, &ione, sol, &ione );
    tvar.real=-dpar[6];
    //r = r - \alpha Ap
    zaxpy(&nd, &tvar, tmp+nd, &ione, tmp+2*nd, &ione );
    //Calculate <r_k+1, r_k+1) 
    zdotc(&tvar, &nd, tmp+2*nd, &ione, tmp+2*nd, &ione );
    dpar[5]=tvar.real;
    //printf(" norm :%f\n",dpar[5]);
    //printf("%f %f\n",tvar.real,tvar.imag);
    ipar[2]=4;
    // Check if <rk,rk> is already small
    *rci=2;
    if (ipar[3]==ipar[4])
      *rci=0;
    return;
  }
  else if (ipar[2]==4){
    //beta = <r_k+1, r_k+1>/<r_k, r_k>
    dpar[7]=dpar[5]/dpar[4];
    tvar.real=dpar[7];
    //p_k+1= r_k+1 + \beta_k p_k
    //FILE *fp;
    //fp=fopen("tempP_prev.dat","w");
    //int i;
    //for(i=0;i<nd;i++)
    //  fprintf(fp,"%f %f\t %f %f\n",tmp[i].real,tmp[i].imag,tmp[2*nd+i].real,tmp[2*nd+i].imag);
    //fclose(fp);
    //printf("saved\n");
    zdscal(&nd, &dpar[7], tmp, &ione);
    //fp=fopen("tempP_scaled.dat","w");
    //for(i=0;i<nd;i++)
    //  fprintf(fp,"%f %f\t %f %f\n",tmp[i].real,tmp[i].imag,tmp[2*nd+i].real,tmp[2*nd+i].imag);
    //fclose(fp);
    //printf("beta=%f\n",dpar[7]);
    zaxpy(&nd, &cone, tmp+2*nd, &ione, tmp, &ione );
    ipar[3]++;
    dpar[4]=dpar[5];
    *rci=1;
    ipar[2]=3;
    return;
  }



}

