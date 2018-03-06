#include "global.h"
#include <stdio.h>

//void generate_fields(){
//  rand_norm_vec(p_mom,nf,0.0,onebyroot2);
//  rand_norm_cvec(R,nf,0.0,onebyroot2);
//  MTx(R,phiu,1);
//}
//
//void measure(){
//  int i;
//  int j;
//  int nbr;
//  double corr_conf=0;
//  for(i=0;i<M;i++){
//    for(j=0;j<nsites;j++){
//      nbr=(j+1)%lx;
//      corr_conf+=2*(Xu[i*nsites+j]*Ru[i*nsites+nbr]+Xu[i*nsites+nbr]*Ru[i*nsites+j]);
//    }
//  }
//  corr_fn+=corr_conf/(double)(1.0*MCS*M*nsites);
//  corr_fn_err+=(corr_conf*corr_conf)/(double)(1.0*MCS*M*M*nsites*nsites);
//}
void calc_deriv(double *deriv, dcomplex *phi){
  int ione=1;
  int block;
  int i;
  int sigma=1;
  dcomplex c1,c2;
  cgsteps+=zconj_grad(X,phi);
  Mx(X,aux3_nf);
  Mprimeax(X,aux2_nf);

  for(block=1;block<M;block++){
    for(i=0;i<nsites;i+=2){
      c1=aux2_nf[i+block*nsites];
      c2=aux3_nf[i+block*nsites];
      deriv[i+(block-1)*nsites]=-2.0*(c1.real*c2.real+c1.imag*c2.imag);
      c1=aux2_nf[(i+1)%lx+block*nsites];
      c2=aux3_nf[(i+1)%lx+block*nsites];
      deriv[i+(block-1)*nsites]+=-2.0*(c1.real*c2.real+c1.imag*c2.imag);
    }

  }
  for(i=0;i<nsites;i+=2){
      c1=aux2_nf[i];
      c2=aux3_nf[i];
      deriv[i+(M-1)*nsites]=2.0*(c1.real*c2.real+c1.imag*c2.imag);
      c1=aux2_nf[(i+1)%lx];
      c2=aux3_nf[(i+1)%lx];
      deriv[i+(M-1)*nsites]+=2.0*(c1.real*c2.real+c1.imag*c2.imag);
  }


  Mprimebx(X,aux2_nf);

  for(block=1;block<M;block++){
    for(i=1;i<nsites;i+=2){
      c1=aux2_nf[i+block*nsites];
      c2=aux3_nf[i+block*nsites];
      deriv[i+(block-1)*nsites]=-2.0*(c1.real*c2.real+c1.imag*c2.imag);
      c1=aux2_nf[(i+1)%lx+block*nsites];
      c2=aux3_nf[(i+1)%lx+block*nsites];
      deriv[i+(block-1)*nsites]+=-2.0*(c1.real*c2.real+c1.imag*c2.imag);
    }

  }
  for(i=1;i<nsites;i+=2){
      c1=aux2_nf[i];
      c2=aux3_nf[i];
      deriv[i+(M-1)*nsites]=2.0*(c1.real*c2.real+c1.imag*c2.imag);
      c1=aux2_nf[(i+1)%lx];
      c2=aux3_nf[(i+1)%lx];
      deriv[i+(M-1)*nsites]+=2.0*(c1.real*c2.real+c1.imag*c2.imag);
  }



}
//Calculate -dV/dx

void calc_force(double *deriv){
  int ione=1;
  int block;
  int i;
  dcomplex c1,c2;
  cgsteps+=zconj_grad(X,phi);
  cg_ctr++;
  Mx(X,aux3_nf);
  Mprimeax(X,aux2_nf);

  for(block=1;block<M;block++){
    for(i=0;i<nsites;i+=2){
      c1=aux2_nf[i+block*nsites];
      c2=aux3_nf[i+block*nsites];
      deriv[i+(block-1)*nsites]=-2.0*(c1.real*c2.real+c1.imag*c2.imag);
      c1=aux2_nf[(i+1)%lx+block*nsites];
      c2=aux3_nf[(i+1)%lx+block*nsites];
      deriv[i+(block-1)*nsites]+=-2.0*(c1.real*c2.real+c1.imag*c2.imag);
    }

  }
  for(i=0;i<nsites;i+=2){
      c1=aux2_nf[i];
      c2=aux3_nf[i];
      deriv[i+(M-1)*nsites]=2.0*(c1.real*c2.real+c1.imag*c2.imag);
      c1=aux2_nf[(i+1)%lx];
      c2=aux3_nf[(i+1)%lx];
      deriv[i+(M-1)*nsites]+=2.0*(c1.real*c2.real+c1.imag*c2.imag);
  }


  Mprimebx(X,aux2_nf);

  for(block=1;block<M;block++){
    for(i=1;i<nsites;i+=2){
      c1=aux2_nf[i+block*nsites];
      c2=aux3_nf[i+block*nsites];
      deriv[i+(block-1)*nsites]=-2.0*(c1.real*c2.real+c1.imag*c2.imag);
      c1=aux2_nf[(i+1)%lx+block*nsites];
      c2=aux3_nf[(i+1)%lx+block*nsites];
      deriv[i+(block-1)*nsites]+=-2.0*(c1.real*c2.real+c1.imag*c2.imag);
    }

  }
  for(i=1;i<nsites;i+=2){
      c1=aux2_nf[i];
      c2=aux3_nf[i];
      deriv[i+(M-1)*nsites]=2.0*(c1.real*c2.real+c1.imag*c2.imag);
      c1=aux2_nf[(i+1)%lx];
      c2=aux3_nf[(i+1)%lx];
      deriv[i+(M-1)*nsites]+=2.0*(c1.real*c2.real+c1.imag*c2.imag);
  }


  for(i=0;i<nf;i++){
    deriv[i]= -dtau*0.25*x_hs[i]-deriv[i];
  }
}

double calc_energy(double *p,double *x){
  double energy=0;
  int i;
  dcomplex c1;
  //free parts
  energy+=dtau*ddot(&nf, x, &ione, x, &ione) +ddot(&nf, p, &ione, p, &ione);
  for(i=0;i<nf;i++)
    aux1_nf[i].real=aux1_nf[i].imag=0;
  zconj_grad(X,phi)
  zdotc(&c1,&nf, phi, &ione, X, &ione);
  energy+=c1.real;
  return energy;

}
//
//
void hamiltonian_evolution(int ifmeasure){
  double t=0;
  int i;
  double e_old,e_new;
  double test;
  //make backups in case metropolis fails;
  dcopy(&nf,p_mom,&ione,pcopy,&ione);
  dcopy(&nf,x_hs,&ione,xcopy,&ione);

  //init momentum : p(t) --> p(t+1/2)
  for(i=0;i<nf;i++)
    X[i].real=X[i].imag=0.0;
  calc_force(dVdx);
  e_old=calc_energy(p_mom,x_hs);

  //measure

  //if(ifmeasure)
  //  measure();

  daxpy(&nf, &dtby2, dVdx, &ione, p_mom, &ione );
  //X =O**(-1) Phi
  for(i=0;i<nf;i++)
    X[i].real=X[i].imag=0.0;
  //Leapfrog
  i=0;
  while(i<=tsteps){
    daxpy(&nf, &twodt, p_mom, &ione, x_hs, &ione);
    //if(i>1){
    //  drotm(&nf, Xu, &ione, Xu_old, &ione, givens_param);
    //  drotm(&nf, Xd, &ione, Xd_old, &ione, givens_param);
    //}
    calc_force(dVdx);
    daxpy(&nf, &dt, dVdx, &ione, p_mom, &ione);
    //if(i==0){
    //  dcopy(&nf,Xu,&ione,Xu_old,&ione);
    //  dcopy(&nf,Xd,&ione,Xd_old,&ione);
    //}
    i++;
  }
  daxpy(&nf, &dt, p_mom, &ione, x_hs, &ione);
  e_new=calc_energy(p_mom,x_hs);

  if((genrand_real2() > exp(e_old-e_new))||(e_new!=e_new)){
    dcopy(&nf,pcopy,&ione,p_mom,&ione);
    dcopy(&nf,xcopy,&ione,x_hs,&ione);
  }
  else{
    accept++;
  }

}








