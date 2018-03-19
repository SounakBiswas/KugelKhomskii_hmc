#include "global.h"
#include <stdio.h>
#include <math.h>
#include "mt19937ar.h"

void rand_norm_vec(double *vec,int dim,double mean,double var);
void make_sparse();
void init_sparse();
//simple one d chain
void make_bonds(){
  int i;
  for (i=0;i<nsites; i++){
    bond2site[i][0]=i;
    bond2site[i][1]=(i+2)%lx;
  }
  for (i=0;i<nsites; i++){
    site2bond[i][0]=i;
    site2bond[i][1]=(lx + i - 2)%lx;
  }
}

int initialize() {
  lx=LX;
  call=0;
  recal=RECAL;
  nsites=NSITES;
  lambdaprime=LAMBDAPRIME;
  twonsites=2*nsites;
  nnbonds=NNBONDS;
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

  init_genrand(26);

  make_bonds();

  Proj=(double *)malloc(nsites*sizeof(double));
  lambda=(int *)malloc((nnbonds+1)*sizeof(int));
  XA=(dcomplex *)malloc(nf*sizeof(dcomplex));
  XB=(dcomplex *)malloc(nf*sizeof(dcomplex));
  XC=(dcomplex *)malloc(nf*sizeof(dcomplex));
  XA_old=(dcomplex*)malloc(nf*sizeof(dcomplex));
  XB_old=(dcomplex*)malloc(nf*sizeof(dcomplex));
  XC_old=(dcomplex*)malloc(nf*sizeof(dcomplex));

  Gr1=(dcomplex*)malloc(nsites*nsites*sizeof(dcomplex));
  Gr2=(dcomplex*)malloc(nsites*nsites*sizeof(dcomplex));

  phiA=(dcomplex *)malloc(nf*sizeof(dcomplex));
  phiB=(dcomplex *)malloc(nf*sizeof(dcomplex));
  phiC=(dcomplex *)malloc(nf*sizeof(dcomplex));

  RA=(dcomplex *)malloc(nf*sizeof(dcomplex));
  RB=(dcomplex *)malloc(nf*sizeof(dcomplex));
  RC=(dcomplex *)malloc(nf*sizeof(dcomplex));


  pcopy=(double *)malloc(nf*sizeof(double));
  xcopy=(double *)malloc(nf*sizeof(double));
  L_pcg=(double *)malloc(nf*sizeof(double));
  D_pcg=(double *)malloc(nf*sizeof(double));
  x_hs=(double *)malloc(nf*sizeof(double));
  aux_pcg=(dcomplex *)malloc(nf*sizeof(dcomplex));
  aux_mtmx_nf=(dcomplex *)malloc(nf*sizeof(dcomplex));
  aux1_nf=(dcomplex *)malloc(nf*sizeof(dcomplex));
  aux2_nf=(dcomplex *)malloc(nf*sizeof(dcomplex));
  aux3_nf=(dcomplex *)malloc(nf*sizeof(dcomplex));
  dVdx=(double *)malloc(nf*sizeof(double));
  aux_ns=(dcomplex*)malloc(nsites*sizeof(dcomplex));
  p_mom=(double *)malloc(nf*sizeof(double));
  rowIndex_kxa=(int*)malloc((nsites+1)*sizeof(int));
  rowIndex_kxb=(int*)malloc((nsites+1)*sizeof(int));
  rowIndex_basis=(int*)malloc((nf+1)*sizeof(int));
  cols_kxa=(int*)malloc(2*nsites*sizeof(int));
  cols_kxb=(int*)malloc(2*nsites*sizeof(int));
  cols_basis=(int*)malloc(2*nf*sizeof(int));

  acsr_kxa=(dcomplex *)calloc(nf*2,sizeof(dcomplex));
  acsr_kxb=(dcomplex *)calloc(nf*2,sizeof(dcomplex));
  acsr_kxap=(dcomplex *)calloc(nf*2,sizeof(dcomplex));
  acsr_kxbp=(dcomplex *)calloc(nf*2,sizeof(dcomplex));
  acsr_basis=(dcomplex *)calloc(nf*2,sizeof(dcomplex));

  //rand_norm_vec(x_hs,nf,0.0,sqrt(2.0/dtau));
  //init_sparse();
  int i;
  for(i=0;i<nnbonds;i++)
    lambda[i]=1;

}
void update_Projs(){
  int i;
  for(i=0;i<nsites;i++)
    Proj[i]=lambda[site2bond[i][0]]*lambda[site2bond[i][1]];
  Proj[0]*=lambdaprime;
  Proj[1]*=lambdaprime;

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
  for(i=0;i<2*nf;i++){
    acsr_kxa[i].imag=acsr_kxa[i].real=0;
    acsr_kxb[i].imag=acsr_kxb[i].real=0;
  }
}
void init_basis(){
  int nbr;
  int link;
  int j;
  int i;
  int ctr;
  ctr=0;
  for(j=0; j<M; j++){
    for(i=0;i<nsites;i++){
      rowIndex_basis[j*nsites+i]=ctr;
      if(i%2==0){
        cols_basis[ctr]=j*nsites+i;
        acsr_basis[ctr].real=1.0/root2;
        acsr_basis[ctr].imag=0;
        ctr++;
        cols_basis[ctr]=j*nsites+(i+1)%lx;
        acsr_basis[ctr].real=1.0/root2;
        acsr_basis[ctr].imag=0;
        ctr++;
      }
      else {
        cols_basis[ctr]=j*nsites+i;
        acsr_basis[ctr].imag=-1.0/root2;
        acsr_basis[ctr].real=0;
        ctr++;
        cols_basis[ctr]=j*nsites+(lx+i-1)%lx;
        acsr_basis[ctr].imag=1.0/root2;
        acsr_basis[ctr].real=0;
        ctr++;
      }

    }
    //getchar();
  }
  rowIndex_basis[nf]=ctr;
}
//populate the hs-field dependent KE matrix
void update_sparse(){
  int nbr;
  int link;
  int j;
  int i;
  int ctr;
  for(j=0; j<M; j++){
    ctr=0;
    //printf("\nblock=%d\n",j);
    for(i=0;i<nsites;i++){
      nbr=(i%2==0)?(i+1)%lx:(lx+i-1)%lx;
      link=(i%2==0)?i:nbr;
      acsr_kxa[j*twonsites+ctr].real=cosh(dtau*x_hs[j*nsites+link]/2.0);
      acsr_kxap[j*twonsites+ctr].real=0.5*dtau*sinh(dtau*x_hs[j*nsites+link]/2.0);

      nbr=(i%2==1)?(i+1)%lx:(lx+i-1)%lx;
      link=(i%2==1)?i:nbr;
      acsr_kxb[j*twonsites+ctr].real=cosh(dtau*x_hs[j*nsites+link]/2.0);
      acsr_kxbp[j*twonsites+ctr].real=0.5*dtau*sinh(dtau*x_hs[j*nsites+link]/2.0);
      ctr++;

      nbr=(i%2==0)?(i+1)%lx:(lx+i-1)%lx;
      link=(i%2==0)?i:nbr;
      if(i%2==0){
        acsr_kxa[j*twonsites+ctr].imag=-sinh(dtau*x_hs[j*nsites+link]/2.0);
        acsr_kxap[j*twonsites+ctr].imag=-dtau*0.5*cosh(dtau*x_hs[j*nsites+link]/2.0);
      }
      else{
        acsr_kxa[j*twonsites+ctr].imag=+sinh(dtau*x_hs[j*nsites+link]/2.0);
        acsr_kxap[j*twonsites+ctr].imag=+0.5*dtau*cosh(dtau*x_hs[j*nsites+link]/2.0);
      }

      nbr=(i%2==1)?(i+1)%lx:(lx+i-1)%lx;
      link=(i%2==1)?i:nbr;
      if(i%2==1){
        acsr_kxb[j*twonsites+ctr].imag=-sinh(dtau*x_hs[j*nsites+link]/2.0);
        acsr_kxbp[j*twonsites+ctr].imag=-dtau*0.5*cosh(dtau*x_hs[j*nsites+link]/2.0);
      }
      else{
        acsr_kxb[j*twonsites+ctr].imag=+sinh(dtau*x_hs[j*nsites+link]/2.0);
        acsr_kxbp[j*twonsites+ctr].imag=+dtau*0.5*cosh(dtau*x_hs[j*nsites+link]/2.0);
      }

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

  free(phiA);
  free(phiB);
  free(phiC);

  free(pcopy);
  free(xcopy);

  free(RA);
  free(RB);
  free(RC);

  free(XA);
  free(XB);
  free(XC);
  free(XA_old);
  free(XB_old);
  free(XC_old);

  free(x_hs);
  free(p_mom);

  free(acsr_kxa);
  free(acsr_kxb);

  free(rowIndex_kxa);
  free(rowIndex_kxb);
  free(cols_kxa);
  free(cols_kxb);

}
