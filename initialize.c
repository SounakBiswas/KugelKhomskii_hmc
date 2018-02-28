#include "global.h"
#include <stdio.h>
#include <math.h>
#include "mt19937ar.h"

void rand_norm_vec(double *vec,int dim,double mean,double var);
void make_sparse();
void init_sparse();

int initialize() {
  lx=LX;
  call=0;
  nsites=NSITES;
  beta=BETA;
  dtau=DTAU;
  M=(int)(beta/dtau);
  nf=nsites*M;
  root2=sqrt(2);
  onebyroot2=1.0/root2;
  root2U=sqrt(2*U);
  dt=T_FL;
  dtby2=dt/2.0;
  twodt=2.0*dt;
  ione=1;
  accept=0;
  tsteps=(int)(1/dt);
  alpha_pcg=1.05;
  cgsteps=0;
  cg_ctr=0;


  corr_fn=0.0;

  givens_param[0]=-1.0;

  givens_param[1]=2.0;
  givens_param[2]=1.0;
  givens_param[3]=-1.0;
  givens_param[4]=0.0;

  init_genrand(25);

  phi=(dcomplex *)malloc(nf*sizeof(double));
  X=(dcomplex *)malloc(nf*sizeof(double));
  X_old=(dcomplex*)malloc(nf*sizeof(double));
  pcopy=(double *)malloc(nf*sizeof(double));
  xcopy=(double *)malloc(nf*sizeof(double));
  L_pcg=(double *)malloc(nf*sizeof(double));
  D_pcg=(double *)malloc(nf*sizeof(double));
  x_hs=(double *)malloc(nf*sizeof(double));
  aux_pcg=(double *)malloc(nf*sizeof(double));
  aux_mtmx_nf=(dcomplex *)malloc(nf*sizeof(dcomplex));
  aux1_nf=(dcomplex *)malloc(nf*sizeof(dcomplex));
  aux2_nf=(dcomplex *)malloc(nf*sizeof(dcomplex));
  aux3_nf=(dcomplex *)malloc(nf*sizeof(dcomplex));
  dVdx=(double *)malloc(nf*sizeof(double));
  R=(dcomplex *)malloc(nf*sizeof(dcomplex));
  aux_ns=(dcomplex*)malloc(nsites*sizeof(dcomplex));
  p_mom=(double *)malloc(nf*sizeof(double));
  rowIndex_kxa=malloc((nsites+1)*sizeof(int));
  rowIndex_kxb=malloc((nsites+1)*sizeof(int));
  cols_kxa=malloc(2*nsites*sizeof(int));
  cols_kxb=malloc(2*nsites*sizeof(int));

  acsr_kxa=(dcomplex *)malloc(nf*2*sizeof(dcomplex));
  acsr_kxb=(dcomplex *)malloc(nf*2*sizeof(dcomplex));

  rand_norm_vec(x_hs,nf,0.0,sqrt(2.0/dtau));
  //init_sparse();

}
void init_sparse(){
  int i,j,k;
  int nnz=2*nsites;
  int ctr;
  ctr=0;
  int nbr;
  for(i=0;i<nsites;i++){
    rowIndex_kxa[i]=ctr;
    rowIndex_kxb[i]=ctr;
    cols_kxa[ctr]=i;
    cols_kxb[ctr]=i;
    ctr++;

    nbr=(i%2==0)?(i+1)%lx:(lx+i-1)%lx;
    cols_kxa[ctr]=nbr;

    nbr=(i%2==1)?(i+1)%lx:(lx+i-1)%lx;
    cols_kxb[ctr]=nbr;

    ctr++;

  }
  rowIndex_kxa[nsites]=ctr;
  rowIndex_kxb[nsites]=ctr;
}
//populate the hs-field dependent KE matrix
void make_sparse(){
  int ctr=0;
  int nbr;
  int link;
  int j;
  int i;
  for(j=0; j<M; j++){
    for(i=0;i<nsites;i++){
      link=(i%2==0)?i:nbr;
      acsr_kxa[M*nsites+ctr].real=cos(dtau*x_hs[link]/2.0);
      ctr++;
      link=(i%2==1)?i:nbr;
      acsr_kxb[M*nsites+ctr].real=cos(dtau*x_hs[link]/2.0);
      ctr++;

      nbr=(i%2==0)?(i+1)%lx:(lx+i-1)%lx;
      link=(i%2==0)?i:nbr;
      if(i%2==0)
        acsr_kxa[M*nsites+ctr].imag=-sinh(dtau*x_hs[link]/2.0);
      else
        acsr_kxa[M*nsites+ctr].imag=sinh(dtau*x_hs[link]/2.0);

      nbr=(i%2==1)?(i+1)%lx:(lx+i-1)%lx;
      link=(i%2==1)?i:nbr;
      if(i%2==0)
        acsr_kxb[M*nsites+ctr].imag=-sinh(dtau*x_hs[link]/2.0);
      else
        acsr_kxb[M*nsites+ctr].imag=sinh(dtau*x_hs[link]/2.0);

      ctr++;
    }
  }
}
void free_all(){
  int i;
  free(aux1_nf);
  free(aux2_nf);
  free(aux3_nf);
  free(aux_ns);
  free(aux_mtmx_nf);

  free(aux_pcg);
  free(L_pcg);
  free(D_pcg);

  free(phi);
  free(pcopy);
  free(xcopy);
  free(R);
  free(x_hs);
  free(p_mom);

  free(acsr_kxa);
  free(acsr_kxb);

  free(rowIndex_kxa);
  free(rowIndex_kxb);
  free(cols_kxa);
  free(cols_kxb);

}
