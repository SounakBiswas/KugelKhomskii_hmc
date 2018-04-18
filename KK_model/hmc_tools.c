#include "global.h"
#include <stdio.h>

void generate_fields(){
  update_Projs();
  update_sparse();
  rand_norm_vec(p_mom,nf,0.0,onebyroot2);
  rand_norm_cvec(RA,nf,0.0,onebyroot2);
  rand_norm_cvec(RB,nf,0.0,onebyroot2);
  rand_norm_cvec(RC,nf,0.0,onebyroot2);
  MDx(RA,phiA);
  MDx(RB,phiB);
  MDx(RC,phiC);
}
//
void measure(){
  int i;
  int j;
  int nbr;
  double corr_conf=0;
  for(i=0;i<M;i++){
    for(j=0;j<nsites;j++){
      nbr=(j+1)%lx;
      corr_conf+=2*(XA[j*nsites+i].real*RA[j*nsites+i].real+XA[j*nsites+i].imag*RA[j*nsites+i].imag);
      //corr_conf+=2*(XB[j*nsites+i].real*RB[j*nsites+i].real+XB[j*nsites+i].imag*RB[j*nsites+i].imag);
    }
  }
  corr_fn+=corr_conf/(double)(1.0*MCS*M*nsites);
  corr_fn_err+=(corr_conf*corr_conf)/(double)(1.0*MCS*M*M*nsites*nsites);
}
void update_phis(){
  cgsteps+=zconj_grad(XA,phiA);
  cgsteps+=zconj_grad(XB,phiB);
  cgsteps+=zconj_grad(XC,phiC);
  cgsteps+=3;

}

//add -d/dx (phi**(dag) (M**D M)**(-1) phi) to deriv
//Use X to store (M**D M)**(-1) phi
void add_dSphidx(double *deriv, dcomplex *phi, dcomplex *X){
  int ione=1;
  int block;
  int i;
  int sigma=1;
  //printf("entered dSphidx \n");
  dcomplex c1,c2;
  //cgsteps+=zconj_grad(X,phi);
  //cg_ctr++;
  //printf("CGdone \n");
  //Mx(X,aux3_nf);
  //Mprimeax(X,aux2_nf);
  Mx(X,aux1_nf);

  //apply projectors.
  /* ***************************** */
  for(i=0;i<0;i++){
    aux1_nf[i].real*=Proj[i];
    aux1_nf[i].real*=Proj[i];
  
  }
  /* ****************************** */

  P2x(aux1_nf,aux3_nf,acsr_kxb,cols_kxb,rowIndex_kxb);

  P1x(X,aux2_nf,acsr_kxap,cols_kxa,rowIndex_kxa);

  for(block=1;block<M;block++){
    for(i=0;i<nsites;i+=2){
      c1=aux2_nf[i+block*nsites];
      c2=aux3_nf[i+block*nsites];
      deriv[i+(block-1)*nsites]+=-2.0*(c1.real*c2.real+c1.imag*c2.imag);
      c1=aux2_nf[(i+1)%lx+block*nsites];
      c2=aux3_nf[(i+1)%lx+block*nsites];
      deriv[i+(block-1)*nsites]+=-2.0*(c1.real*c2.real+c1.imag*c2.imag);
    }

  }
  for(i=0;i<nsites;i+=2){
      c1=aux2_nf[i];
      c2=aux3_nf[i];
      deriv[i+(M-1)*nsites]+=2.0*(c1.real*c2.real+c1.imag*c2.imag);
      c1=aux2_nf[(i+1)%lx];
      c2=aux3_nf[(i+1)%lx];
      deriv[i+(M-1)*nsites]+=2.0*(c1.real*c2.real+c1.imag*c2.imag);
  }


  //Mprimebx(X,aux2_nf);
  Mx(X,aux3_nf);
  P1x(X,aux1_nf,acsr_kxa,cols_kxa,rowIndex_kxa);
  P2x(aux1_nf,aux2_nf,acsr_kxbp,cols_kxb,rowIndex_kxb);
  //apply projectors.
  /* ***************************** */
  for(i=0;i<0;i++){
    aux2_nf[i].real*=Proj[i];
    aux2_nf[i].real*=Proj[i];
  
  }


  for(block=1;block<M;block++){
    for(i=1;i<nsites;i+=2){
      c1=aux2_nf[i+block*nsites];
      c2=aux3_nf[i+block*nsites];
      deriv[i+(block-1)*nsites]+=-2.0*(c1.real*c2.real+c1.imag*c2.imag);
      c1=aux2_nf[(i+1)%lx+block*nsites];
      c2=aux3_nf[(i+1)%lx+block*nsites];
      deriv[i+(block-1)*nsites]+=-2.0*(c1.real*c2.real+c1.imag*c2.imag);
    }

  }
  for(i=1;i<nsites;i+=2){
      c1=aux2_nf[i];
      c2=aux3_nf[i];
      deriv[i+(M-1)*nsites]+=2.0*(c1.real*c2.real+c1.imag*c2.imag);
      c1=aux2_nf[(i+1)%lx];
      c2=aux3_nf[(i+1)%lx];
      deriv[i+(M-1)*nsites]+=2.0*(c1.real*c2.real+c1.imag*c2.imag);
  }



}
//Calculate -dV/dx

void calc_force(double *deriv){
  int i;
  for(i=0;i<nf;i++){
    deriv[i]= -dtau*0.25*x_hs[i];
  }
  add_dSphidx(deriv,phiA,XA);
  add_dSphidx(deriv,phiB,XB);
  add_dSphidx(deriv,phiC,XC);


}

double calc_energy(double *p,double *x){
  double energy=0;
  double temp;
  int i;
  dcomplex c1;
  //free parts
   temp=0.125*dtau*ddot(&nf, x, &ione, x, &ione) +ddot(&nf, p, &ione, p, &ione);
  energy+=temp;
//  printf("free part=%f\t",temp);
  for(i=0;i<nf;i++)
    aux1_nf[i].real=aux1_nf[i].imag=0;
  //zconj_grad(X,phi);
  zdotc(&c1,&nf, phiA, &ione, XA, &ione);
  //printf("A part=%f\t",c1.real);
  energy+=c1.real;

  zdotc(&c1,&nf, phiB, &ione, XB, &ione);
  //printf("B part=%f\n",c1.real);
  energy+=c1.real;

  zdotc(&c1,&nf, phiC, &ione, XC, &ione);
  energy+=c1.real;
  return energy;

}
void hamiltonian_evolution(int ifmeasure){
  double t=0;
  int i;
  double e_old,e_new;
  double test;
  update_Projs();
  update_sparse();
  //make backups in case metropolis fails;
  dcopy(&nf,p_mom,&ione,pcopy,&ione);
  dcopy(&nf,x_hs,&ione,xcopy,&ione);

  //init momentum : p(t) --> p(t+1/2)
  for(i=0;i<nf;i++){
    XA[i].real=XA[i].imag=0.0;
    XB[i].real=XB[i].imag=0.0;
    XC[i].real=XC[i].imag=0.0;
  }
  update_phis();
  calc_force(dVdx);
  e_old=calc_energy(p_mom,x_hs);
  printf("energy old:%f\t",e_old);

  //measure

  if(ifmeasure)
    measure();

  daxpy(&nf, &dtby2, dVdx, &ione, p_mom, &ione );
  //X =O**(-1) Phi
  for(i=0;i<nf;i++){
    XA[i].real=XA[i].imag=0.0;
    XB[i].real=XB[i].imag=0.0;
    XC[i].real=XC[i].imag=0.0;
  }
  //Leapfrog
  i=0;
  while(i<=tsteps){
    daxpy(&nf, &twodt, p_mom, &ione, x_hs, &ione);
    update_sparse();
    if(i>1){
      zdrotm( XA,  XA_old,  givens_param);
      zdrotm( XB,  XB_old,  givens_param);
      zdrotm( XC,  XC_old,  givens_param);
    }
    update_phis();
    calc_force(dVdx);
    daxpy(&nf, &dt, dVdx, &ione, p_mom, &ione);
    if(i==0){
      zcopy(&nf,XA,&ione,XA_old,&ione);
      zcopy(&nf,XB,&ione,XB_old,&ione);
      zcopy(&nf,XC,&ione,XC_old,&ione);
    }
    i++;
  }
  daxpy(&nf, &dt, p_mom, &ione, x_hs, &ione);
  e_new=calc_energy(p_mom,x_hs);
  printf("energy new:%f\n",e_new);

  if((genrand_real2() > exp(e_old-e_new))||(e_new!=e_new)){
    dcopy(&nf,pcopy,&ione,p_mom,&ione);
    dcopy(&nf,xcopy,&ione,x_hs,&ione);
    //printf("not \n");
  }
  else{
    //printf("accept \n");
    accept++;
  }

}








