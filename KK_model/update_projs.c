#include "global.h"
#include <stdio.h>
void det_update(){
  int i,j;
  dcomplex alpha,beta=0;
  alpha.real=1;
  alpha.im=beta.im=beta.real=0;
  //init to identity
  for (i=0;i<nsites; i++){
    for (j=0;j<nsites; j++){
      if(i==j)
        Gr1[i+j*nsites].real=1;
      else
        Gr1[i+j*nsites].real=0;
      Gr1[i+j*nsites].imag=0;
    }
  }
  char trans[1]="N";
  char matdescra="MLNCAA";

  //no stabilization
  for(i=0;i<M;i++){
    zcsrmm (trans , &nsites , &nsites , &nsites , &alpha , matdescra , acsr_kxa ,  , rowIndex_kxa , rowIndex_kxa+1 , Gr1 , &nsites , &beta , &Gr2 , &nsites );
    zcsrmm (trans , &nsites , &nsites , &nsites , &alpha , matdescra , acsr_kxb ,  , rowIndex_kxb , rowIndex_kxb+1 , Gr2 , &nsites , &beta , &Gr1 , &nsites );
  }
  // Multiply projector matrices, add identity
  for(i=0;i<nsites; i++){
    for (j=0;j<nsites; j++){
      Gr1[i+j*nsites]*=Proj[i];
      if(i==j)
        Gr1[i+j*nsites]+=1;
    }
  }
  invert_matrix(Gr1,nsites);
  int site1,site2;
  double del;
  double ratio1, ratio2;

  for(i=0;i<nsites;i++){
    site1=bond2site[i][0];
    site2=bond2site[i][1];
    del=-2*lambda[i];
    ratio1=1+(1-Gr1[site1+nsites*site1])*del;
    for(j=0;j<nsites;j++)
      for(k=0;k<nsites;k++)
        Gr2[j+k*nsites]=Gr1[j+k*nsites] - Gr1[j+k*nsites]*del*((site1==k) - Gr1[i+k*nsites])/ratio1;

    ratio2=1+(1-Gr2[site2+nsites*site2])*del;
    if(genrand_real2()<(ratio1*ratio2)){
      for(j=0;j<nsites;j++)
        for(k=0;k<nsites;k++)
          Gr1[j+k*nsites]=Gr1[j+k*nsites] - Gr2[j+k*nsites]*del*((site1==k) - Gr2[i+k*nsites])/ratio2;
      lambda[i]*=-1;
    }
  }
}
void cg_update(){
  double del;
  double ratio1, ratio2;
  double eold,enew;
  update_Projs();
  update_sparse();
  update_phis();
  eold=calc_energy(p_mom,x_hs);
  //loop over bonds

  for(i=0;i<nsites;i++){
    lambda[i]*=-1;
    update_Projs();
    update_sparse();
    update_phis();
    enew=calc_energy(p_mom,x_hs);
    if(genrand_real2()>exp(eold-enew)){
      //revert
      lambda[i]*=-1;
    }
    else{
      //update e_old
      eold=enew;
    }

  }





}
