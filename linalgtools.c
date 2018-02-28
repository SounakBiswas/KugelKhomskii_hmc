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
void MTx(dcomplex *x,dcomplex *y){
  int block;
  int ione=1;
  dcomplex minusone;
  dcomplex one;
  minusone.real=-1.0;
  minusone.imag=0.0;
  one.real=1.0;
  one.imag=0.0;
  for(block=0;block<(M-1);block++){
    zcsrsymv("L",&nsites, acsr_kxb+block*nsites,rowIndex_kxb, cols_kxb, x+(block+1)*nsites, aux_ns);
    zcsrsymv("L",&nsites, acsr_kxa+block*nsites,rowIndex_kxa, cols_kxa, aux_ns, y+block*nsites);
    zscal (&nsites, &minusone, y+block*nsites, &ione);
    zaxpy (&nsites, &one, x+block*nsites,&ione, y+block*nsites, &ione);
  }
  zcsrsymv("L",&nsites, acsr_kxb +(M-1)*nsites,rowIndex_kxb, cols_kxb, x, aux_ns);
  zcsrsymv("L",&nsites, acsr_kxa +(M-1)*nsites,rowIndex_kxa, cols_kxa, aux_ns,y+(M-1)*nsites);
  zscal(&nsites, &one, y+(M-1)*nsites, &ione);
  zaxpy(&nsites, &one, x+(M-1)*nsites, &ione, y+(M-1)*nsites, &ione);
}
void Mx(dcomplex *x,dcomplex *y){
  int block;
  int i;
  int ione=1;
  dcomplex minusone;
  dcomplex one;
  minusone.real=-1.0;
  minusone.imag=0.0;
  one.real=1.0;
  one.imag=0.0;
  zcsrsymv("L", &nsites, acsr_kxa+M*nsites, rowIndex_kxa, cols_kxa,x+(M-1)*nsites,aux_ns);
  zcsrsymv("L", &nsites, acsr_kxb+M*nsites, rowIndex_kxb, cols_kxb,aux_ns,y);
  zaxpy(&nsites, &one, x, &ione, y, &ione);
  for(block=1;block<M;block++){
    zcsrsymv ("L",&nsites, acsr_kxa+(block-1)*nsites,rowIndex_kxa, cols_kxa,x+(block-1)*nsites,aux_ns);
    zcsrsymv ("L",&nsites, acsr_kxb+(block-1)*nsites,rowIndex_kxb, cols_kxb, aux_ns, y+block*nsites);
    zscal (&nsites, &minusone, y+block*nsites,&ione);
    zaxpy (&nsites, &one ,x+block*nsites, &ione, y+block*nsites, &ione);
  }
}
//void Mxprime(double*x,double *y,int sigma){
//  int block;
//  int i;
//  int ione=1;
//  double minusone=-1;
//  double one=1;
//  dcsrsymv("L",&nsites, acsr_kxa,rowIndex_kxa, cols_kxa,x+(M-1)*nsites,aux_ns);
//  dcsrsymv("L",&nsites, acsr_kxb,rowIndex_kxb, cols_kxb,aux_ns,y);
//  for(block=1;block<M;block++){
//    dcsrsymv ("L",&nsites, acsr_kxa,rowIndex_kxa, cols_kxa,x+(block-1)*nsites,aux_ns);
//    dcsrsymv ("L",&nsites, acsr_kxb,rowIndex_kxb, cols_kxb,aux_ns,y+block*nsites);
//  }
//}

void MTMx(dcomplex *x,dcomplex *y){
  Mx(x,aux_mtmx_nf);
  MTx(aux_mtmx_nf,y);
}
//calc_force= -dV/dx









