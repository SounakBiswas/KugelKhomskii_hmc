#include "global.h"
#include "math.h"
#include "mkl.h"
#include <assert.h>
#include <stdio.h>

int zconj_grad(dcomplex *sol,dcomplex *rhs){

  int i,j,k;
  int  n=nf;
  int ipar[128];
  double dpar[128];
  dcomplex *tmp;
  double *tmp1;
  tmp=(dcomplex*)malloc(nf*4*sizeof(dcomplex));
  double norm;
  double minusone=-1.0;
  double one=1;
  int rci_request;
  //MKL_INT ione=1;
  double norm_rhs=dznrm2(&n,rhs,&ione);
  //printf("norm_rhs=%f\n",norm_rhs);

  //init solution
  //for(i=0;i<nf;i++)
  //  sol[i]=0.0;

  zcg_init (n, sol, rhs, &rci_request, ipar, dpar, tmp);
  //ipar[7]=0;
  ipar[4]=1000;
  ipar[10]=IFPRECON;
  //if(rci_request==0)
    //printf("initialized\n");

  //dcg_check (&n, sol, rhs, &rci_request, ipar, dpar, tmp);
  //if(rci_request==0)
    //printf("checked\n");

  rci_request=1;
  while(rci_request) {
    zcg (n, sol, rhs, &rci_request, ipar, dpar, tmp);
//    printf("req=%d %d\n",rci_request,ipar[2]);
//    getchar();
    if(rci_request==0)
      break;
    if(rci_request==1)
      MDMx (tmp, tmp+nf) ;
    if(rci_request==3)
      apply_preconditioner (tmp+2*nf, tmp+3*nf) ;

    if(rci_request==2){
      norm=dznrm2(&n,tmp+2*nf,&ione);
      //printf("%f %f %16f\n",norm,norm_rhs,(norm/norm_rhs));
      //getchar();
      if((norm/norm_rhs)< 1.E-6)
        break;
    }
    //if(rci_request==3)
    //  apply_preconditioner(tmp+2*nf,tmp+3*nf,sigma);
  }
  //MKL_INT itercount;
  //dcg_get (&n, sol, rhs, &rci_request, ipar, dpar, tmp,&itercount);
  free(tmp);

  mkl_free_buffers();
  //printf("steps %d\n",ipar[3]);
  return ipar[3];


}
