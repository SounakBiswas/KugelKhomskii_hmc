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
  rand_norm_vec(x_hs,nf,0.0,sqrt(16.0/dtau));
  init_basis();
  init_sparse();
  //printf("wassup \n");
  //FILE*fp;
  //fp=fopen("vecr.dat","r");
  //for(i=0;i<nf;i++)
  //  fscanf(fp,"%lf\n",&phiA[i].real);
  //fclose(fp);
  //fp=fopen("veci.dat","r");
  //for(i=0;i<nf;i++)
  //  fscanf(fp,"%lf\n",&phiA[i].imag);
  //fclose(fp);
  //fp=fopen("hs.dat","r");
  //for(i=0;i<nf;i++)
  //  fscanf(fp,"%lf\n",&x_hs[i]);
  //fclose(fp);
  update_sparse();
  update_Projs();
  for(i=0;i<nsites;i++)
    printf("%f\t",Proj[i]);
  getchar();

  for(i=0;i<WUP;i++){
    if(i%50==0)
    printf("%d\n",i);
    generate_fields();
    //printf("gen fields done\n");
    hamiltonian_evolution(0);
  }
  printf("acceptance rate= %f\n", (double)(accept)/(double)(1.0*(WUP)));
  printf("avg cg iterations= %f\n", (double)(cgsteps)/(double)(1.0*(cg_ctr)));
  cgsteps=0;
  cg_ctr=0;
  accept=0;
  for(i=0;i<MCS;i++){
    if(i%50==0)
    printf("%d\n",i);
    generate_fields();
    //printf("gen fields done\n");
    hamiltonian_evolution(1);
  }
  printf("avg cg iterations= %f\n", (double)(cgsteps)/(double)(1.0*(cg_ctr)));
  printf("acceptance rate= %f\n", (double)(accept)/(double)(1.0*(MCS)));
  printf("corr_fn= : %f %f\n", corr_fn,sqrt((corr_fn_err-corr_fn*corr_fn)/MCS));
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

  //for(i=0;i<nf;i++)
  //  dVdx[i]=0;
  //MDMx(phi,aux1_nf);
  //printf("enter conj grad\n");
  //dcomplex tvar;
  //zdotc(&tvar, &nf, phi,&ione, phi, &ione);
  //double t2=dznrm2(&nf, phi,&ione);
  //printf("%f %f %f\n",t2,tvar.real, tvar.imag);
  //zconj_grad(aux1_nf,phi);
  //add_dSphidx(dVdx,phiA);
  //hamiltonian_evolution(0);
  //
  //fp=fopen("vect.dat","w");
  //for(i=0;i<nf;i++)
  //  fprintf(fp,"%lf\t%lf\n",phi[i].real,phi[i].imag);
  //fclose(fp);
  //fp=fopen("vec3.dat","w");
  //for(i=0;i<nf;i++)
  //  fprintf(fp,"%lf\n",dVdx[i]);
  //fclose(fp);
}

