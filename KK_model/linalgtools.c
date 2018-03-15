#include "global.h"
#include "math.h"
#include "mkl.h"
#include <stdio.h>

/***********************************************************************************************/


/*calling dsyev subroutine to diagonalize matrix x of dimension y. The array xeigenval[] would
  contain the eigenvalues and the array xeigenvec[] would contain the eigenvectors. dsyev implimentation of lapack
  is used for our case*/


//y=M**T x
//thse routines multiply M and MT using the checkerboard decomposition
void MDx(dcomplex *x,dcomplex *y){
  int block;
  int ione=1;
  int i;
  char iftransp[1]={'N'};
  dcomplex minusone;
  dcomplex one;
  minusone.real=-1.0;
  minusone.imag=0.0;
  one.real=1.0;
  one.imag=0.0;
  for(block=0;block<(M-1);block++){
    zcsrgemv(iftransp,&nsites, acsr_kxb+block*twonsites,rowIndex_kxb, cols_kxb, x+(block+1)*nsites, aux_ns);
    zcsrgemv(iftransp,&nsites, acsr_kxa+block*twonsites,rowIndex_kxa, cols_kxa, aux_ns, y+block*nsites);
    zscal (&nsites, &minusone, y+block*nsites, &ione);
    zaxpy (&nsites, &one, x+block*nsites,&ione, y+block*nsites, &ione);
  }
  //using y as a temp buffer first to enable the projector
  for(i=0;i <nsites;i++){
    y[i+(M-1)*nsites].real=Proj[i]*x[i].real;
    y[i+(M-1)*nsites].imag=Proj[i]*x[i].imag;
  }
  zcsrgemv(iftransp,&nsites, acsr_kxb +(M-1)*twonsites,rowIndex_kxb, cols_kxb, y+(M-1)*nsites, aux_ns);
  zcsrgemv(iftransp,&nsites, acsr_kxa +(M-1)*twonsites,rowIndex_kxa, cols_kxa, aux_ns,y+(M-1)*nsites);
  zscal(&nsites, &one, y+(M-1)*nsites, &ione);
  zaxpy(&nsites, &one, x+(M-1)*nsites, &ione, y+(M-1)*nsites, &ione);
}
void Mx(dcomplex *x,dcomplex *y){
  int block;
  int i;
  int ione=1;
  char iftransp[1]={'N'};
  dcomplex minusone;
  dcomplex one;
  minusone.real=-1.0;
  minusone.imag=0.0;
  one.real=1.0;
  one.imag=0.0;
  zcsrgemv(iftransp, &nsites, acsr_kxa+(M-1)*twonsites, rowIndex_kxa, cols_kxa,x+(M-1)*nsites,aux_ns);
  zcsrgemv(iftransp, &nsites, acsr_kxb+(M-1)*twonsites, rowIndex_kxb, cols_kxb,aux_ns,y);
  for(i=0;i <nsites;i++){
    y[i].real*=Proj[i];
    y[i].imag*=Proj[i];
  }
  zaxpy(&nsites, &one, x, &ione, y, &ione);
  for(block=1;block<M;block++){
    zcsrgemv (iftransp,&nsites, acsr_kxa+(block-1)*twonsites,rowIndex_kxa, cols_kxa,x+(block-1)*nsites,aux_ns);
    zcsrgemv (iftransp,&nsites, acsr_kxb+(block-1)*twonsites,rowIndex_kxb, cols_kxb, aux_ns, y+block*nsites);
    zscal (&nsites, &minusone, y+block*nsites,&ione);
    zaxpy (&nsites, &one ,x+block*nsites, &ione, y+block*nsites, &ione);
  }
}
// M- I = P2.P1
void P1x(dcomplex *x,dcomplex *y,dcomplex *acsr, int *cols, int *rowIndex){
  int block;
  int i;
  int ione=1;
  char iftransp[1]={'N'};
  dcomplex minusone;
  dcomplex one;
  minusone.real=-1.0;
  minusone.imag=0.0;
  one.real=1.0;
  one.imag=0.0;
  zcsrgemv(iftransp, &nsites, acsr+(M-1)*twonsites, rowIndex, cols,x+(M-1)*nsites,y);
  for(block=1;block<M;block++){
    zcsrgemv (iftransp,&nsites, acsr+(block-1)*twonsites,rowIndex, cols,x+(block-1)*nsites,y+block*nsites);
  }
}
void P2x(dcomplex *x,dcomplex *y,dcomplex *acsr, int *cols, int *rowIndex){
  int block;
  int i;
  int ione=1;
  char iftransp[1]={'N'};
  dcomplex minusone;
  dcomplex one;
  minusone.real=-1.0;
  minusone.imag=0.0;
  one.real=1.0;
  one.imag=0.0;
  zcsrgemv(iftransp, &nsites, acsr+(M-1)*twonsites, rowIndex, cols,x,y);
  for(block=1;block<M;block++){
    zcsrgemv (iftransp,&nsites, acsr+(block-1)*twonsites,rowIndex, cols,x+block*nsites,y+block*nsites);
  }
}
void MDMx(dcomplex *x,dcomplex *y){
  Mx(x,aux_mtmx_nf);
  MDx(aux_mtmx_nf,y);
}
void zdrotm(dcomplex *X, dcomplex *X_old, double* givens_param){
  int i;
  dcomplex temp;
  for(i=0; i<nf; i++){
    temp=X[i];
    X[i].real= 2*X[i].real - X_old[i].real;
    X[i].imag= 2*X[i].imag - X_old[i].imag;
    X_old[i]=temp;
  }


}
//calc_force= -dV/dx









