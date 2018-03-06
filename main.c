#include <stdio.h>
#include <stdlib.h>
#include "math.h"
#include "global.h"
#include "mkl.h"

void initialize();
void init_fields();
void mkl_free_buffers(void);
void Mx(dcomplex *,dcomplex*);
void make_sparse();

//void main(){
//  int i,j,k;
//  initialize();
//  for(i=0;i<WUP;i++){
//    generate_fields();
//    hamiltonian_evolution(0);
//  }
//  printf("warmup done, avg cg iterations= %f\n",(double)(cgsteps)/(double)(1.0*(cg_ctr)));
//  printf("acceptance rate= %f\n", (double)(accept)/(double)(1.0*(WUP)));
//  for(i=0;i<MCS;i++){
//    generate_fields();
//    hamiltonian_evolution(1);
//  }
//  printf("acceptance rate= %f\n", (double)(accept)/(double)(1.0*(MCS+WUP)));
//  printf("corr_fn= : %f %f\n", corr_fn,sqrt((corr_fn_err-corr_fn*corr_fn)/MCS));
//  printf("avg cg iterations= %f\n", (double)(cgsteps)/(double)(1.0*(cg_ctr)));
//  free_all();
//}
void main(){
  int i,j,k;
  initialize();
  init_sparse();
  dcomplex *dense;
  dense=(dcomplex *)malloc(nsites*nsites*sizeof(dcomplex));
  FILE*fp;
  fp=fopen("hs.dat","r");
  for(i=0;i<nf;i++)
    fscanf(fp,"%lf\n",&x_hs[i]);
  fclose(fp);
  fp=fopen("vecr.dat","r");
  for(i=0;i<nf;i++)
    fscanf(fp,"%lf\n",&phi[i].real);
  fclose(fp);
  fp=fopen("veci.dat","r");
  for(i=0;i<nf;i++)
    fscanf(fp,"%lf\n",&phi[i].imag);
  fclose(fp);
  make_sparse();
  MKL_INT job[5];
  job[0]=1;
  job[1]=0;
  job[2]=0;
  job[3]=2;
  MKL_INT info;
  //mkl_zdnscsr(job,&nsites, &nsites, dense, &nsites, acsr_kxb+(0)*twonsites,cols_kxb,rowIndex_kxb,&info);
  //for(i=0; i<nsites; i++){
  //  for(j=0; j<nsites; j++)
  //    printf("%f+i%f\t",dense[i+j*nsites].real,dense[i+j*nsites].imag);
  //  printf("\n");
  //}
  //getchar();
  //mkl_zdnscsr(job,&nsites, &nsites, dense, &nsites, acsr_kxa+(0)*twonsites,cols_kxa,rowIndex_kxa,&info);
  //for(i=0; i<nsites; i++){
  //  for(j=0; j<nsites; j++)
  //    printf("%f+i%f\t",dense[i+j*nsites].real,dense[i+j*nsites].imag);
  //  printf("\n");
  //}
  //getchar();

  for(i=0;i<nf;i++)
    aux1_nf[i].real=aux1_nf[i].imag=0;
  //MDMx(phi,aux1_nf);
  //printf("enter conj grad\n");
  //dcomplex tvar;
  //zdotc(&tvar, &nf, phi,&ione, phi, &ione);
  //double t2=dznrm2(&nf, phi,&ione);
  //printf("%f %f %f\n",t2,tvar.real, tvar.imag);
  //zconj_grad(aux1_nf,phi);
  calc_deriv(aux1_nf,phi);
  fp=fopen("vect.dat","w");
  for(i=0;i<nf;i++)
    fprintf(fp,"%lf\t%lf\n",phi[i].real,phi[i].imag);
  fclose(fp);
  fp=fopen("vec3.dat","w");
  for(i=0;i<nf;i++)
    fprintf(fp,"%lf\t%lf\n",aux1_nf[i].real,aux1_nf[i].imag);
  fclose(fp);
}

